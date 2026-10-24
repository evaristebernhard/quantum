# DQI Reed–Solomon BM reproduction on Kaggle

This folder reproduces the published Berlekamp–Massey (BM) reversible decoder used in the DQI Reed–Solomon task.

## What runs

- Kaggle kernel: `varistebernhard/rs-bm-reproduction`
- Notebook: `main.ipynb`
- Upstream implementation: [andigu/uncompute](https://github.com/andigu/uncompute)
- Pinned upstream commit: `64b5a13aa249f1360846657392e6c66817b26d0f`
- Evaluator: upstream `tools/eval.py dqi-reed-solomon`
- Expected baseline metrics: 442,397 Toffoli gates, 1,837 peak live qubits, TQ-cost 812,683,289.
- The evaluator checks phase, syndrome preservation, and workspace cleanup under the task's test distribution.

The circuit is an approximate construction under the paper's stated failure budget; it is not a proof of worst-case exact decoding for every possible input. This reproduction is a baseline, not a new optimization.

## CPU-only Kaggle run and result retrieval

The kernel is configured as private with `enable_gpu: "false"`; no `machine_shape` is requested. Internet access remains enabled because the notebook fetches the pinned upstream source and installs its locked dependencies.

GitHub Actions uses the GitHub Environment named `KAGGLE_API_TOKEN` and reads the Environment secret of the same name via `secrets.KAGGLE_API_TOKEN`. The token itself is never written to this repository.

After pushing the notebook, the workflow waits for the Kaggle run and then attempts to download its outputs, even if the notebook reports an execution failure. The notebook writes these files under `/kaggle/working` so Kaggle can expose them as kernel outputs:

- `summary.json` — compact resource metrics and verdict
- `bm-reproduction-report.json` — evaluator report
- `run-manifest.json` — pinned source revision and command
- `evaluator-stdout.txt` and `evaluator-stderr.txt` — logs
- `bm-reproduction-results.zip` — a bundle of result files

The workflow uploads the downloaded files as the `rs-bm-reproduction-results` artifact (retained for 30 days). To retrieve them, open the run under GitHub Actions and download the artifact. If the Kaggle download step fails, the artifact also contains a diagnostic file and workflow metadata.

## Trigger

Pushes to `master` or `main` that change this folder or the workflow trigger the push. You can also start a run manually from GitHub Actions → Run workflow.

## Attribution

The notebook fetches and runs source code from [andigu/uncompute](https://github.com/andigu/uncompute) at the pinned commit above. That upstream repository is licensed Apache-2.0. No upstream implementation files are copied into this repository.
