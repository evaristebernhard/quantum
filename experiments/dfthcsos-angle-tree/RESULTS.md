# Angle-tree screen results

**Status: research screen only; no official physical verification or score was run.**

Run: [GitHub Actions 40165322832](https://github.com/evaristebernhard/quantum/actions/runs/40165322832)  
Artifact: `dfthcsos-angle-tree-screen` (expires 2026-11-01).

## Scope

The screen used the repository's starter DFTHC-SOS factorization, not the separate low-cost R1/B719/C350 research candidate discussed in the literature review. The starter shape is R=65, B=60, C=60, with 3,976 angle rows. Angles were re-encoded to 8-bit signed-gauge words. All results below are weighted projector-distortion proxies, not the official error_1/error_2, SOS gap, normalization, Toffoli count, or competition score.

## Tree-topology search

The script enumerated all one-edge swaps from the current topology, ranked them by the first-order sensitivity proxy, then evaluated deterministic weighted projector distortion for the 12 shortlisted candidates. The workflow repeated this for up to three greedy rounds, retaining only a candidate that improved the deterministic proxy.

| Checkpoint | Weighted projector distortion | Relative change |
|---|---:|---:|
| Reordered starter-tree baseline | 0.0932924193 | — |
| After round 1 | 0.0916564562 | -1.754% |
| After round 2 | 0.0912060909 | -0.491% |
| After round 3 | 0.0904682755 | -0.809% |
| **Total** | **0.0904682755** | **-3.027%** |

Accepted topology changes:

| Round | Removed edge | Added edge |
|---|---|---|
| 1 | (2, 34) | (34, 59) |
| 2 | (6, 38) | (32, 38) |
| 3 | (4, 36) | (14, 36) |

The first-order proxy is only weakly predictive of the deterministic result: the best deterministic candidate ranked 4th, 10th, and 3rd by the proxy in rounds 1–3. That is evidence to keep the proxy as a shortlist generator, not as the acceptance criterion.

## Edge-order/QROM screen

The second stage sampled 128 valid topological orderings. It preserved parent-before-child and sibling order, so only commuting rotations were reordered and their angle-word columns permuted together.

- Active difference payload: 49 bits in every one of the 11 batches before and after.
- Recursive QROM-node Toffoli proxy: 32,626 before and after.
- No ordering improvement was found under this proxy.

This is a useful negative result: for this starter factorization, simple edge reordering does not reduce the active-bit width or the measured recursive QROM proxy. The current word tables are already saturated on those metrics.

## Conclusion and next gate

The tree-topology method is mechanically working and reduces the deterministic projector proxy by about 3%, but it has **not** demonstrated a reduction in official physical error or scored circuit cost. The edge-order-only route is not promising on this starter candidate.

The next meaningful test is to run the same screen on a physically valid low-cost factorization and include actual official error/gap metrics plus built ANGLE-circuit costs for a small Pareto shortlist. Do not promote any NPZ from this screen to a submission without matching circuit files and a passing unchanged verifier.
