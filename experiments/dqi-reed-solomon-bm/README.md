# DQI Reed–Solomon BM reproduction on Kaggle

This folder reproduces the published Berlekamp–Massey (BM) reversible decoder used in the DQI Reed–Solomon task.

## What runs

- Kaggle kernel: `varistebernhard/rs-bm-reproduction`
- Notebook: `main.ipynb`
- Upstream implementation: [andigu/uncompute](https://github.com/andigu/uncompute)
- Pinned upstream commit: `64b5a13aa249f1360846657392e6c66817b26d0f`
- Evaluator: upstream `tools/eval.py dqi-reed-solomon`
- Expected resource metrics at this revision: 442,397 Toffoli gates, 1,837 peak live qubits, TQ-cost 812,683,289.
- The evaluator checks the output phase, unchanged syndrome registers, and workspace cleanup under the task's test distribution.

The circuit is an approximate construction under the paper's stated failure budget; it is not a proof of worst-case exact decoding for every possible input. This reproduction must not be described as a new optimization.

## GitHub Actions / Kaggle setup

1. In GitHub repository Settings → Environments, identify the Environment that contains the secret `KAGGLE_API_TOKEN`.
2. In `.github/workflows/kaggle-push.yml`, replace `your_environment_name` with that exact Environment name.
3. Ensure `KAGGLE_API_TOKEN` is an Environment **secret** (not a committed file or plain repository variable).
4. Merge the workflow and experiment files into `master` (the current default branch) or `main`. Pushes to either branch that touch this experiment trigger `kaggle kernels push`; the push starts a Kaggle run. A manual run is also available via GitHub Actions → Run workflow.
5. Kaggle must permit GPU and Internet access for the account. The kernel is private and requests `NvidiaTeslaT4`, which is the Kaggle identifier for T4 ×2.

## GPU caveat

The T4 ×2 accelerator is configured and the notebook reports whether it is present. The upstream circuit builder and resource evaluator are CPU-oriented, however, and do not dispatch their gate simulation to CUDA. The GPU is therefore allocated but is not expected to speed up this reproduction. The useful output is the reproducible circuit cost and correctness verdict, not GPU utilization.

## Attribution

The notebook fetches and runs source code from [andigu/uncompute](https://github.com/andigu/uncompute) at the pinned commit above. That upstream repository is licensed Apache-2.0. No upstream implementation files are copied into this repository.
