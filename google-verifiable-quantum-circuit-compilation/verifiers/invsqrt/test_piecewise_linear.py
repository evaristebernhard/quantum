from __future__ import annotations

import numpy as np

import piecewise_linear_circuit as candidate


def test_coefficient_table_is_complete_and_fits_registers() -> None:
    assert len(candidate.COEFFICIENTS) == 64
    assert all(len(row) == 2 for row in candidate.COEFFICIENTS)
    assert all(0 < value < (1 << candidate.COEF_BITS)
               for row in candidate.COEFFICIENTS for value in row)


def test_exhaustive_benchmark_domain_relative_error() -> None:
    n_p = 10
    max_x = 3 * (2**n_p - 2)**2
    xs = np.arange(1, max_x + 1, dtype=np.int64)
    errors = candidate.relative_error(n_p, xs)
    worst_index = int(np.argmax(errors))
    print(
        f"max relative error={errors[worst_index]:.12e}, "
        f"worst x={int(xs[worst_index])}, target={2**-24:.12e}"
    )
    assert np.isfinite(errors).all()
    assert errors.max() < 2**-24
