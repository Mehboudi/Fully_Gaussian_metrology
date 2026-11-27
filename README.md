# Fully_Gaussian_metrology

MATLAB code for the Monte Carlo simulation and maximum likelihood estimation (MLE) in the squeezing-estimation example of

> J. Navarro, S. Morelli, M. Sanz, and M. Mehboudi, *On super additivity of Fisher information in fully Gaussian metrology*, DOI: [10.1088/2058-9565/aec629](https://doi.org/10.1088/2058-9565/aec629).
Codes published in Zenodo: https://doi.org/10.5281/zenodo.24592598

The script reproduces **Fig. 2** (main text) and **Fig. 4** (Appendix E) of the paper.

## Contents

- `Simulation_of_outcomes_and_MLE.m`: the main MATLAB script. Running it opens two figure windows:
  1. **"MLE Analysis: Convergence and Scaling"** → **Fig. 2** of the paper (2×2 panels):
     - (a) top-left: one MLE trajectory for $m=2$, converging to the true value of $r$ as the number of repetitions $\nu$ grows
     - (b) top-right: Fisher information per copy vs. the number of copies $m$, comparing the global Gaussian strategy, the local Gaussian strategy, and the QFI
     - (c) bottom-left: MSE of the MLE vs. $\nu$ with the CRB (global), the QCRB and the CRB (local), for $m=2$
     - (d) bottom-right: the same for $m=5$
  2. **"MLE Analysis: MSE Breakdown"** → **Fig. 4** of the paper (Appendix E): the squared bias (solid) and the variance (dashed) of the MLE vs. $\nu$ for $m=2$ (blue) and $m=5$ (red). The shaded bands are 95% bootstrap percentile confidence intervals for the squared bias.
- `LICENSE`: MIT license.

## Requirements

- MATLAB (base MATLAB only: `randn`, `randi`, `fminsearch`, `accumarray`, `fill`, `waitbar`, ...). No toolboxes are needed.

## Usage

1. Open MATLAB and set the current folder to the repository root.
2. If needed, edit the parameters at the top of `Simulation_of_outcomes_and_MLE.m`:

   | Parameter | Default | Meaning |
   |-----------|---------|---------|
   | `alpha0`  | `sqrt(2)` | Amplitude $\alpha$ of the coherent probe, $d=[\alpha,\alpha]^T$ |
   | `nu`      | `100`     | Number of repetitions $\nu$ (sequential measurement rounds) per trajectory |
   | `N_sims`  | `100000`  | Number of independent Monte Carlo trajectories |
   | `true_r`  | `0.1`     | True squeezing parameter $r$ used to generate the data |
   | `B_boot`  | `1000`    | Number of bootstrap resamples for the squared-bias confidence bands (set further down, in the bootstrap section) |

3. Run the script from the MATLAB command window:
```matlab
   run Simulation_of_outcomes_and_MLE.m
```
4. To save the figures, add these lines at the end of the script, or run them afterwards:
```matlab
   exportgraphics(figure(1), 'Fig2_MLE_convergence.png', 'Resolution', 300);
   exportgraphics(figure(2), 'Fig4_bias_variance.png',  'Resolution', 300);
```
   (`savefig` works too if you want `.fig` files.)

## What the script does

**Measurement model.** The $m$ copies are measured jointly with the global strategy from the paper. A linear optical network moves the whole displacement onto the $m$-th mode. That mode is homodyned in the quadrature that is optimal for the displacement. The other $m-1$ modes are homodyned in the quadrature that is optimal for the covariance matrix. For each of the $\nu$ rounds the simulated outcomes are therefore

- $x^{(j)}_k \sim \mathcal{N}\!\left(0,\ e^{2r}/2\right)$ for $j = 1,\dots,m-1$,
- $x^{(m)}_k \sim \mathcal{N}\!\left(\sqrt{m/2}\,\alpha\,(e^{r}-e^{-r}),\ \cosh(2r)/2\right)$.

**Estimator.** For every trajectory and every $\nu' = 1,\dots,\nu$, the MLE is found with `fminsearch` by maximising the log-likelihood (Eq. (E2) of the paper). The estimate from the previous step is used as the starting point:

$$
\log p = -\sum_{k=1}^{\nu'}\left[\frac{\sum_{j=1}^{m-1}\big(x^{(j)}_k\big)^2}{e^{2r}} + \frac{\big(\sqrt{m/2}\,\alpha\,(e^{r}-e^{-r}) - x^{(m)}_k\big)^2}{\cosh(2r)}\right] - \nu'(m-1)\,r - \frac{\nu'}{2}\log\cosh(2r) + \text{const.}
$$

**Error analysis.** From the $N$ trajectories the script computes, at each $\nu$:

- the mean squared error $\mathrm{MSE}(\hat r) = \langle(\hat r - r)^2\rangle$, which is plotted in Fig. 2(c, d);
- the bias $b(\hat r) = \langle \hat r\rangle - r$ and the variance $\mathrm{Var}(\hat r)$, which are plotted in Fig. 4. Recall that $\mathrm{MSE} = \mathrm{Var} + b^2$.
- a 95% bootstrap confidence band for $b^2$. The script resamples the $N$ whole trajectories with replacement $B$ times, recomputes $b^2$ for every $\nu$, and takes the 2.5% and 97.5% percentiles.

**Theoretical bounds** (Fisher information per round of $m$ copies, evaluated around $r_0 = 0$):

- Global Gaussian strategy: $\mathcal{F}^C_{\rm Global} = 4m\alpha^2 + 2(m-1)$
- Optimal local Gaussian strategy: $\mathcal{F}^C_{\rm Local} = 4m\alpha^2$
- Quantum Fisher information: $\mathcal{F}^Q = m(4\alpha^2 + 2)$

The bounds in panels (c, d) are $1/(\nu\,\mathcal{F})$.

## Code structure

- `run_simulation(m, nu, N_sims, alpha0, true_r)` simulates the outcomes and runs the sequential MLE. It returns:
  - `all_estimates` (`nu × N_sims`): the estimate of $r$ after each round, for every trajectory
  - `squared_errors` (`nu × N_sims`): $(\hat r - r)^2$
  - `avg_error` (`nu × 1`): the MSE vs. $\nu$
  - `bias` (`nu × 1`): $\langle \hat r\rangle - r$ vs. $\nu$
  - `est_var` (`nu × 1`): the variance of $\hat r$ over trajectories vs. $\nu$
- `neg_log_likelihood(r, X_m, S1, nu_curr, m, alpha)` is the negative log-likelihood minimised by `fminsearch`.
- `bootstrap_sqbias(all_estimates, true_r, B)` returns the lower and upper 95% bootstrap percentile bounds on the squared bias for each $\nu$.
- `shade_band(x, lo, hi, col)` draws the shaded confidence band on the log–log axes.

## Runtime and memory

- The default `N_sims = 100000` gives the statistics used in the paper, but the run is long. The script calls `fminsearch` `nu × N_sims` times for each of $m = 2$ and $m = 5$. For quick tests, use `N_sims = 1000` or fewer. With very few trajectories the bias estimate and its bootstrap band will be noisy.
- The script stores every trajectory (`nu × N_sims` arrays for the estimates and squared errors, for both values of $m$), which takes a few hundred MB with the defaults. The bootstrap step needs these full arrays, so you can't replace them with running averages if you want Fig. 4.

## Citation

If you use this code, please cite the paper:

> J. Navarro, S. Morelli, M. Sanz, and M. Mehboudi, *On super additivity of Fisher information in fully Gaussian metrology*, DOI: [10.1088/2058-9565/aec629](https://doi.org/10.1088/2058-9565/aec629).

## License

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
