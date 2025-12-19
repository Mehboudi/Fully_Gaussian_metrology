# Fully_Gaussian_metrology

Simulations and data analysis for Maximum Likelihood Estimation (MLE) in squeezing estimation.

This repository contains MATLAB code that reproduces Figure 2 from our preprint, entitled "\\maketitle ".

## Contents

- `Simulation_of_outcomes_and_MLE.m` — main MATLAB script. Running this script generates the 2x2 panel figure (Figure 2 in the preprint) showing:
  - Panel 1 (top-left): single MLE trajectory for m=2
  - Panel 2 (top-right): Fisher information (per sample) vs. number of modes m
  - Panel 3 (bottom-left): MSE convergence and theoretical bounds for m=2
  - Panel 4 (bottom-right): MSE convergence and theoretical bounds for m=5

## Requirements

- MATLAB (tested with recent versions; the code uses only base MATLAB functions such as `randn`, `fminsearch`, `figure`, and `waitbar`)

No additional toolboxes are strictly required to run the provided script.

## Usage

1. Open MATLAB and set the current folder to the repository root (or add it to the MATLAB path).
2. Edit the top-of-file parameter block in `Simulation_of_outcomes_and_MLE.m` if you want to change simulation settings:
   - `alpha0` — tunable parameter alpha
   - `nu` — total number of samples per simulation (number of sequential measurements)
   - `N_sims` — number of independent Monte Carlo simulations (default is large; reduce for quick tests)
   - `true_r` — the ground-truth parameter value used to generate synthetic data
3. Run the script:
   - In the MATLAB command window: `run Simulation_of_outcomes_and_MLE.m`
4. The script displays a 2x2 figure. To save the figure (recommended), add one of the following lines at the end of the script:
   - Save as PNG: `exportgraphics(gcf, 'Figure2.png', 'Resolution', 300);`
   - Save as MATLAB figure: `savefig('Figure2.fig');`

## Notes on performance and suggestions

- The default `N_sims = 100000` is chosen for good Monte Carlo statistics but will increase runtime. For development and quick checks, try `N_sims = 1000` or smaller.
- Memory: the script stores all estimates as an array sized `nu x N_sims`. For very large `N_sims` this can grow large; consider computing running averages or saving only statistics instead of every trajectory if memory becomes an issue.

## Reproducing Figure 2

The file `Simulation_of_outcomes_and_MLE.m` contains everything required to reproduce the figure. Key parameters that control the figure are at the top of the file; make sure you set them to the values used in the paper. After running, save the produced figure using `exportgraphics` or `saveas` for inclusion in the manuscript.

## Script overview (quick)

- `run_simulation(m, nu, N_sims, alpha0, true_r)` — runs N_sims Monte Carlo experiments for a measurement configuration with `m` modes and returns:
  - `all_estimates` (nu-by-N_sims): estimated parameter r at each sample
  - `squared_errors` (nu-by-N_sims)
  - `avg_error` (nu-by-1): MSE vs number of samples (averaged over simulations)
- `neg_log_likelihood(r, X_m, S1, nu_curr, m, alpha)` — negative log-likelihood used by `fminsearch` to find the MLE at each step

## Citation / Authors

Please cite the corresponding preprint / paper when using these results. (DOI:)

## License

(Add license here — e.g., MIT, BSD — or state "All rights reserved" if not open source.)
