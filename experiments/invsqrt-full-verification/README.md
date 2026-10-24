# Full invsqrt verification on Kaggle CPU

This experiment takes the piecewise-linear inverse-square-root circuit candidate and runs the repository's official KickMix verifier on Kaggle with **GPU disabled**.

## Full workload

- Benchmark: Google Verifiable Quantum Circuit Compilation, Benchmark I
- Input domain: every x from 1 through 3,133,452
- Shots: 64 per x (200,540,928 total simulated shots)
- Checks: phase tolerance, input-register restoration, scratch-qubit cleanup, and official resource score
- Reference comparison: circuit qubit/gate counts and the core-cost proxy are saved separately; the proxy is not substituted for the official score

## Automation

The `.github/workflows/invsqrt-kaggle.yml` workflow pushes `main.ipynb` to Kaggle, waits for the CPU notebook to finish, downloads `/kaggle/working` outputs, and uploads them as the `invsqrt-full-verification-results` GitHub Actions artifact.

The job uses the GitHub Environment named `KAGGLE_API_TOKEN` and its `KAGGLE_API_TOKEN` secret. The token is never stored in the repository.

The full verifier is intentionally not reduced to a small input range. If the Kaggle runtime fails or times out, the workflow still attempts to retrieve logs and partial outputs.
