"""Experimental 32-segment piecewise-linear invsqrt circuit.

This is a side-by-side candidate for the Google Verifiable Quantum Circuit
Compilation Benchmark I. The reference circuit is left unchanged.

Approximation:
  * M in [1, 2): 32 intervals of width 1/32.
  * M in [2, 4): 32 intervals of width 1/16.
  * y0 = C0 - T(C1 * delta), with 17 fractional bits.
  * One Newton correction and the reference circuit's remaining arithmetic.

The coefficient lookup is implemented as an algebraic-normal-form (ANF) Boolean
network. Its monomials are shared across all 34 coefficient output bits.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT_DIR = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(PROJECT_DIR / "kickmix"))

import nr_reference_circuit as nr  # noqa: E402

km = nr.km

# Rows 0..31 cover M in [1, 2); rows 32..63 cover M in [2, 4).
# Each pair is (C0, C1), with coefficients represented in units of 2^-17.
COEFFICIENTS: tuple[tuple[int, int], ...] = (
    (131064, 64053), (129063, 61186), (127151, 58527), (125322, 56101),
    (123570, 53817), (121888, 51643), (120274, 49653), (118722, 47769),
    (117229, 46016), (115791, 44387), (114404, 42799), (113066, 41331),
    (111774, 39945), (110526, 38675), (109318, 37433), (108149, 36257),
    (107016, 35104), (105919, 34084), (104854, 33038), (103821, 32080),
    (102818, 31170), (101844, 30336), (100896, 29466), (99975, 28693),
    (99078, 27914), (98206, 27229), (97355, 26494), (96527, 25866),
    (95719, 25216), (94931, 24578), (94163, 24030), (93412, 23425),
    (92677, 22652), (91262, 21638), (89910, 20701), (88617, 19844),
    (87377, 19011), (86189, 18279), (85047, 17557), (83950, 16906),
    (82894, 16284), (81877, 15694), (80897, 15156), (79950, 14605),
    (79037, 14139), (78154, 13677), (77300, 13241), (76473, 12818),
    (75672, 12415), (74896, 12045), (74144, 11703), (73413, 11353),
    (72704, 11037), (72015, 10730), (71345, 10434), (70693, 10139),
    (70059, 9870), (69442, 9624), (68841, 9382), (68255, 9144),
    (67684, 8921), (67127, 8706), (66583, 8484), (66053, 8304),
)

COEF_BITS = 17
DELTA_BITS = 14
SQUARE_BITS = 28
PRODUCT_BITS = 28
RESULT_BITS = 29
INDEX_BITS = 5
ADDRESS_BITS = 6
MONOMIAL_MASKS = tuple(mask for mask in range(1, 1 << ADDRESS_BITS)
                       if mask.bit_count() >= 2)
MONOMIAL_INDEX = {mask: i for i, mask in enumerate(MONOMIAL_MASKS)}


def _mobius_transform(values: list[int]) -> tuple[int, ...]:
    """Return the ANF coefficients of a Boolean function on six input bits."""
    result = list(values)
    for bit in range(ADDRESS_BITS):
        for mask in range(1 << ADDRESS_BITS):
            if mask & (1 << bit):
                result[mask] ^= result[mask ^ (1 << bit)]
    return tuple(result)


ANF_TABLES = tuple(
    tuple(
        _mobius_transform([(row[column] >> output_bit) & 1 for row in COEFFICIENTS])
        for output_bit in range(COEF_BITS)
    )
    for column in range(2)
)


def fixed_point_model(n_p: int, xs: np.ndarray) -> np.ndarray:
    """Reproduce the candidate circuit's integer arithmetic for every x in xs.

    The normalization, square, M multiplication, Newton correction, output
    scaling, and row-wise truncation follow nr_reference_circuit.py. Only the
    initial cubic Horner stage is replaced by a piecewise-linear seed.
    """
    if n_p < 7:
        raise ValueError("The current 14-bit delta layout requires n_p >= 7.")
    width = 2 * n_p + 2
    xs = np.asarray(xs, dtype=np.int64)
    if not ((1 <= xs) & (xs < 2**width)).all():
        raise ValueError("x must lie in [1, 2^width).")

    def bit(v: np.ndarray, i: int) -> np.ndarray:
        return (v >> i) & 1

    def product(rows: list[tuple[np.ndarray, int]], operand: np.ndarray,
                out_bits: int, operand_bits: int, limit: int) -> np.ndarray:
        acc = np.zeros_like(xs)
        for control, weight in rows:
            shift = out_bits - operand_bits + weight
            if shift >= 0:
                acc += control * (operand << shift)
            else:
                acc += control * (operand >> -shift)
        if (acc >= 2**limit).any():
            raise OverflowError("Intermediate register overflow.")
        return acc

    # Normalize x exactly as the reference circuit does.
    s = (width - (np.floor(np.log2(xs)).astype(np.int64) + 1)) // 2
    x_norm = xs << (2 * s)
    branch = bit(x_norm, width - 1).astype(np.int64)
    low_index = (x_norm >> (2 * n_p - 5)) & 31
    high_index = (x_norm >> (2 * n_p - 4)) & 31
    index = np.where(branch == 0, low_index, high_index).astype(np.int64)
    address = index + 32 * branch
    table = np.asarray(COEFFICIENTS, dtype=np.int64)
    c0 = table[address, 0]
    c1 = table[address, 1]

    # delta is M minus the interval's left endpoint, truncated to 14 bits.
    m_q14 = x_norm >> (2 * n_p - DELTA_BITS)
    base_q14 = np.where(branch == 0,
                        (1 << DELTA_BITS) + (index << 9),
                        (2 << DELTA_BITS) + (index << 10))
    delta_q14 = m_q14 - base_q14
    # Match the circuit's column-truncated Horner rows exactly: each delta
    # bit contributes (c1 >> m) independently; do not multiply then truncate once.
    horner_sub = np.zeros_like(xs)
    for m in range(5, DELTA_BITS + 1):
        delta_bit = (delta_q14 >> (DELTA_BITS - m)) & 1
        horner_sub += delta_bit * (c1 >> m)
    y0 = c0 - horner_sub
    if ((y0 <= 0) | (y0 >= 2**COEF_BITS)).any():
        raise OverflowError("Initial approximation does not fit its register.")

    # Match nr._square's column-truncated partial products. A single integer
    # multiplication followed by >> 6 can carry discarded low-column terms into
    # the retained result, while the actual circuit discards each term separately.
    y_bits = {k: bit(y0, COEF_BITS - k) for k in range(1, COEF_BITS + 1)}
    square = np.zeros_like(xs)
    for k in range(1, COEF_BITS + 1):
        if 2 * k <= SQUARE_BITS:
            square += y_bits[k] << (SQUARE_BITS - 2 * k)
        for ell in range(k + 1, COEF_BITS + 1):
            if k + ell - 1 <= SQUARE_BITS:
                square += (y_bits[k] & y_bits[ell]) << (SQUARE_BITS - (k + ell - 1))

    # T(M * square), preserving the reference's per-row truncation semantics.
    m_rows = [(bit(x_norm, i), i - 2 * n_p) for i in range(width)]
    z = product(m_rows, square, PRODUCT_BITS, SQUARE_BITS, PRODUCT_BITS + 2)
    correction = (3 << PRODUCT_BITS) - z

    # T(y0 * correction), again using the reference's row-wise truncation.
    y_bits = {k: bit(y0, COEF_BITS - k) for k in range(1, COEF_BITS + 1)}
    p = product([(y_bits[k], -k) for k in range(1, COEF_BITS + 1)],
                correction, RESULT_BITS, PRODUCT_BITS, RESULT_BITS + 1)
    return p << s


def relative_error(n_p: int, xs: np.ndarray) -> np.ndarray:
    """Classical-model relative error against 1/sqrt(x)."""
    xs = np.asarray(xs, dtype=np.int64)
    output = fixed_point_model(n_p, xs)
    approx = output / 2.0 ** (RESULT_BITS + n_p + 1)
    exact = 1.0 / np.sqrt(xs.astype(np.float64))
    return np.abs(approx - exact) / exact


def _select_index(builder: km.CircuitBuilder, x: km.array,
                  index: km.array, scratch: km.array) -> None:
    """XOR the five-bit segment index into index; the operation is self-inverse.

    For M in [1, 2), index bits come from x[2*n_p-1:2*n_p-5].
    For M in [2, 4), they come from x[2*n_p:2*n_p-4].
    The branch bit selects between the two bit slices.
    """
    width = len(x)
    branch = x[width - 1]
    diff, and_tmp = scratch[0], scratch[1]
    for j in range(INDEX_BITS):
        low = x[width - 7 + j]
        high = x[width - 6 + j]
        builder.cx(low, index[j])
        builder.cx(low, diff)
        builder.cx(high, diff)
        nr._and(builder, False, branch, diff, and_tmp)
        builder.cx(and_tmp, index[j])
        nr._and(builder, True, branch, diff, and_tmp)
        builder.cx(low, diff)
        builder.cx(high, diff)


def _load_coefficients(builder: km.CircuitBuilder, x: km.array,
                       index: km.array, coefs: list[km.array],
                       scratch: km.array) -> None:
    """XOR-load C0 and C1 using a shared six-input ANF lookup network.

    Address bits are index[0:5] (least significant first), followed by the
    branch bit. All nonlinear monomials are computed once and shared by every
    output bit, then uncomputed in reverse dependency order.
    """
    address = [index[i] for i in range(INDEX_BITS)] + [x[len(x) - 1]]
    monomials = {mask: scratch[MONOMIAL_INDEX[mask]] for mask in MONOMIAL_MASKS}

    for mask in MONOMIAL_MASKS:
        lowest = mask & -mask
        bit_index = lowest.bit_length() - 1
        rest = mask ^ lowest
        if rest.bit_count() == 1:
            source = address[rest.bit_length() - 1]
        else:
            source = monomials[rest]
        nr._and(builder, False, source, address[bit_index], monomials[mask])

    for column, register in enumerate(coefs):
        for output_bit in range(COEF_BITS):
            for mask, enabled in enumerate(ANF_TABLES[column][output_bit]):
                if not enabled:
                    continue
                if mask == 0:
                    builder.x(register[output_bit])
                elif mask.bit_count() == 1:
                    builder.cx(address[mask.bit_length() - 1], register[output_bit])
                else:
                    builder.cx(monomials[mask], register[output_bit])

    for mask in reversed(MONOMIAL_MASKS):
        lowest = mask & -mask
        bit_index = lowest.bit_length() - 1
        rest = mask ^ lowest
        if rest.bit_count() == 1:
            source = address[rest.bit_length() - 1]
        else:
            source = monomials[rest]
        nr._and(builder, True, source, address[bit_index], monomials[mask])


def _linear_horner(builder: km.CircuitBuilder, inverse: bool, x: km.array,
                   selectors: km.array, coefs: list[km.array],
                   scratch: km.array) -> None:
    """Compute c0 -= T(c1 * delta), with delta in [0, 1/16).

    The m=5 delta bit is present only in the high branch; m=6..14 use the same
    x bits in both branches. selector[5] is a temporary AND and is cleared
    before this stage returns.
    """
    width = len(x)
    branch = x[width - 1]
    delta5 = selectors[5]
    nr._and(builder, False, branch, x[width - 7], delta5)
    c0, c1 = coefs
    for m in range(5, DELTA_BITS + 1):
        control = delta5 if m == 5 else x[width - 2 - m]
        terms = [(c1[i], i - m) for i in range(m, len(c1))]
        nr._add_row(builder, inverse, control, terms, c0, scratch, subtract=True)
    nr._and(builder, True, branch, x[width - 7], delta5)


def _compute_y1(builder: km.CircuitBuilder, inverse: bool, prec: nr.Precision,
                x: km.array, flags: km.array, selectors: km.array,
                coefs: list[km.array], square: km.array, z: km.array,
                y1: km.array, scratch: km.array) -> None:
    y0 = coefs[0]
    p = y1[:len(y1) - len(flags)]
    with builder.mark("linear_uncompute_y1" if inverse else "linear_compute_y1"):
        if inverse:
            nr._denormalize(builder, True, flags, y1)
            nr._multiply_y0(builder, True, y0, z, p, scratch)
            nr._correct(builder, True, z)
            nr._multiply_m(builder, True, x, square, z, scratch)
            nr._square(builder, True, y0, square, scratch)
            _linear_horner(builder, True, x, selectors, coefs, scratch)
            _load_coefficients(builder, x, selectors[:INDEX_BITS], coefs, scratch)
            _select_index(builder, x, selectors[:INDEX_BITS], scratch)
            nr._normalize(builder, True, x, flags)
        else:
            nr._normalize(builder, False, x, flags)
            _select_index(builder, x, selectors[:INDEX_BITS], scratch)
            _load_coefficients(builder, x, selectors[:INDEX_BITS], coefs, scratch)
            _linear_horner(builder, False, x, selectors, coefs, scratch)
            nr._square(builder, False, y0, square, scratch)
            nr._multiply_m(builder, False, x, square, z, scratch)
            nr._correct(builder, False, z)
            nr._multiply_y0(builder, False, y0, z, p, scratch)
            nr._denormalize(builder, False, flags, y1)


def build_circuit(n_p: int = 10, *, compute_only: bool = False) -> km.Circuit:
    """Build the piecewise-linear candidate circuit."""
    if n_p < 7:
        raise ValueError(f"n_p must be at least 7, got {n_p}.")
    width = 2 * n_p + 2
    prec = nr.Precision()
    out_bits = prec.output_bits(n_p)
    builder = km.CircuitBuilder()

    def alloc(length: int, name: str) -> km.array:
        if compute_only:
            return builder.create_quantum_register(length, name=name)
        return builder.alloc_qubits(length)

    x = builder.create_quantum_register(width, name="x")
    y1 = alloc(out_bits, "y1")
    flags = alloc(n_p, "flags")
    selectors = alloc(INDEX_BITS + 1, "selectors")
    coefs = [alloc(COEF_BITS, "c0"), alloc(COEF_BITS, "c1")]
    square = alloc(SQUARE_BITS, "square")
    z = alloc(PRODUCT_BITS + 2, "z")
    scratch = alloc(max(SQUARE_BITS, PRODUCT_BITS + 2, RESULT_BITS + 1,
                        len(MONOMIAL_MASKS)), "scratch")

    _compute_y1(builder, False, prec, x, flags, selectors, coefs,
                square, z, y1, scratch)
    if not compute_only:
        nr._phase_by_bits(builder, y1)
        _compute_y1(builder, True, prec, x, flags, selectors, coefs,
                    square, z, y1, scratch)
        for register in (y1, flags, selectors, *coefs, square, z, scratch):
            builder.free(register)
    return builder.finish_circuit()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", help="Output .kmb or .kmx circuit file.")
    parser.add_argument("--n-p", type=int, default=10)
    parser.add_argument("--format", choices=("kmb", "kmx"), default="kmb")
    parser.add_argument("--skip-model", action="store_true",
                        help="Skip the exhaustive classical-model error scan.")
    args = parser.parse_args()

    max_x = 3 * (2**args.n_p - 2)**2
    if not args.skip_model:
        errors = relative_error(args.n_p, np.arange(1, max_x + 1, dtype=np.int64))
        worst_x = int(np.argmax(errors)) + 1
        print(f"Classical model: max relative error {errors.max():.12e} at x={worst_x}; "
              f"tolerance={2**-24:.12e}.")
        if errors.max() >= 2**-24:
            raise SystemExit("Classical fixed-point model exceeds the target tolerance.")

    circuit = build_circuit(args.n_p)
    circuit.to_file(args.output, format=args.format)
    stats = circuit.operation_statistics()
    print(f"Wrote {args.output} ({args.format}): {circuit.num_qubits} qubits, "
          f"{stats.num_toffoli_like:,} Toffoli-like operations, {len(circuit):,} operations.")


if __name__ == "__main__":
    main()
