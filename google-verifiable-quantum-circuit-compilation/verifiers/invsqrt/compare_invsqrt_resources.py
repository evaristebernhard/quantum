from __future__ import annotations

import argparse
import json

import nr_reference_circuit as reference
import piecewise_linear_circuit as candidate


def describe(circuit):
    stats = circuit.operation_statistics()
    return {
        "qubits": int(circuit.num_qubits),
        "operations": int(len(circuit)),
        "toffoli_like": int(stats.num_toffoli_like),
        "t_like": int(stats.num_t_like),
        "z_pow_level_4_or_higher": int(
            stats.num_z_pow(min_clifford_hierarchy_level=4)
        ),
        "core_cost_proxy_Q_times_Toffoli_plus_T_over_4": (
            int(circuit.num_qubits)
            * (float(stats.num_toffoli_like) + float(stats.num_t_like) / 4.0)
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-p", type=int, default=10)
    args = parser.parse_args()
    baseline = reference.build_nr_circuit(args.n_p)
    proposed = candidate.build_circuit(args.n_p)
    report = {
        "n_p": args.n_p,
        "baseline_newton_raphson": describe(baseline),
        "candidate_piecewise_linear": describe(proposed),
        "note": (
            "The proxy excludes synthesized non-Clifford Z rotations and is not "
            "the official verifier score. Run the verifier for an authoritative score."
        ),
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
