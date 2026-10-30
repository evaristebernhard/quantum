# Public frontier reproduction, 2026-10-30

This experiment imports hash-pinned public circuits, not a new independently designed compiler. Full provenance and component-specific licensing are in CREDITS.md and LICENSE.md. In particular, the DQI lineage has an upstream license caveat; a compiled artifact does not constitute delivery of a generator for winner eligibility.

## Artifact source

Public source: https://www.kaggle.com/datasets/heliosli/qcc-circuits-live, dataset version 23. QCC_RELEASE_MANIFEST.json pins all 18 submission file hashes. package-manifest.json pins the ZIP submitted to the competition. Large downloaded/generated binaries and the local official verifier snapshot are intentionally not committed.

## Environment and verifier

Existing Conda environment kaggle, Python 3.11, official PyPI kickmix 1.0.0. run-manifest.json records the repository revision used for the unmodified official verifier snapshot. No pre-existing BM edits were changed.

## Quick validation only

- invsqrt: 62 logical qubits, cost 101002.53550915231, ratio 20.5616049463. Eight workers, 12.0 seconds, one shot per x in 1..3133440, followed by a 256-shot boundary check over the 12 highest valid inputs. This is NOT the official 64 measurements per input. Worst sampled phase error/tolerance is 0.99932377 at x=161; synthesis margin is accounted for by the verifier.
- DQI: 889 logical qubits, 4096 shots, 0 value and 0 phase errors, four workers, 30.8 seconds. Mean Toffoli count 433785.2070. Default confidence -128 retained; its low-sample repetition factor is 2.9847756. Quick cost 1.023269301239e12, quick ratio 23.8837064.
- Chemistry: 430 logical qubits, 16266 Toffolis, one shot per address (1220 shots), 9.9 seconds, cost 337993966.8795, ratio 9.92798664. Lambda_eff 48.32364940; both errors remain below 20 mHa.

The quick-sample geometric mean is 16.95664676. A HYPOTHETICAL 102400-shot DQI run with zero errors and unchanged average resource counts would give DQI ratio 59.8420801 and portfolio projection 23.0307612. This assumption is not verification, and neither figure is a Kaggle leaderboard result. Actual scoreboard scoring must be read from the submission result.

The first invsqrt command in run_quick.py was rejected by local verifier parameter validation because its shot count was not divisible by batch size 256. inv_quick_aligned.py corrected the batching and completed successfully. Historical failure reports are retained for audit.

## Files and reproduction

run_quick.py and inv_quick_aligned.py contain the exact quick-check commands. quick-summary-final.json contains raw quick results and the explicitly hypothetical DQI projection. candidate/ contains only the official ZIP members. Do not publish a claim of full verification based on these quick reports.
