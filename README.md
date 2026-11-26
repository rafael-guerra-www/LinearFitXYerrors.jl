# LinearFitXYerrors.jl

<img src="https://user-images.githubusercontent.com/20739393/132685450-1a34351f-ad02-49ee-9d57-b498437e356a.png" width="320" />

This small Julia package, based on York et al. (2004), performs 1D linear fitting of experimental data with uncertainties in both X and Y:

            Linear fit:             Y = a + b*X                             [1]
            
            Errors:                 X ± σX;  Y ± σY                         [2]

            Errors' correlation:    r =  = cov(σX, σY) / (σX * σY)          [3]

where:
- `X` and `Y` are input data vectors with length ≥ 3
- Optional standard deviation errors `σX` and `σY` are vectors or scalars
- Optional `r` is the correlation between the `σX` and `σY` errors\
           `r` can be a vector or scalar

For the `σX` and `σY` errors (error ellipses) a bivariate Gaussian distribution is assumed.\
If no errors are provided, or if only `σX` or `σY` are provided, then the results are equivalent to those obtained using the [LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl) package.

The package computes:
- The intercept `a`, the slope `b` and their uncertainties `σa` and `σb`
- `σa95` and `σb95`: 95%-confidence interval uncertainties corrected by two-tailed t-Student distribution, e.g.: `b ± σb95 = b ± t(0.975,N-2)*σb`
- Goodness of fit `S` (reduced Χ² test): the underlying quantity has Χ² distribution with N-2 degrees of freedom\
  `S ~ 1`: fit consistent with errors, `S > 1`: poor fit, `S >> 1`: errors underestimated, `S < 1`: errors overestimated
- Pearson's correlation coefficient `ρ` of the observed data (standard,
  descriptive convention, comparable with other software)
- Optional display of results with error ellipses and the 95% confidence
  band of the fitted line (`plot(st)` or `isplot=true`)

## Plotting

Since v0.2.0, plotting is optional and backend-friendly: `Plots.jl` is a
*weak* dependency loaded automatically as a package extension when available
(requires Julia ≥ 1.9). Without `Plots`, the package loads fast and everything
except plotting works.

```julia
using Plots                    # optional; activates the extension
using LinearFitXYerrors
st = linearfitxy(X, Y; σX, σY, r)
plot(st)                       # error ellipses + 95% CI band of the line
```

The figure shows the fit line with its pointwise 95% confidence band
(`level = 0.95` by default), per-point error bars (or covariance-error
ellipses when both X and Y have correlated errors), and the fit summary in
the title. The legacy `isplot=true` keyword still works via the extension.

##
## Installation
```julia
julia> ] add LinearFitXYerrors
julia> using LinearFitXYerrors
```
##
## Usage
```julia
# Optional: plotting (loads the Plots extension)
# ] add Plots

# The input data and regression results are returned in the fields of the `st` structure (::stfitxy):

st = linearfitxy(X, Y)    # no errors in X and Y, no plot displayed

st = linearfitxy(X, Y; σX, σY, isplot=true)    # X-Y errors not correlated (r=0); plot with ratio=1

st = linearfitxy(X, Y; σX, σY, r=0, isplot=true, ratio=\:auto)  # X-Y errors not correlated (r=0); plot with ratio=1

# Or plot afterwards (Plots.jl loaded):
# st = linearfitxy(X, Y; σX, σY, r)
# plot(st)
```

## Notes:
- The objective for this first package was to learn how to publish a Julia package via Github while implementing York's technique.
- Since v0.2.0 the console summary of the fit is opt-in: pass `verbose=true` to `linearfitxy` to print it (previously always printed).
- The 95% confidence band of the fitted line is drawn for all input cases (data with or without errors).
- The package author is not a statistician and the topics of "errors in variables" and "confidence intervals" are beyond his expertise.
- While the results are consistent with the references provided, the apparent exception of Amen (2012) is now understood (see the note on Amen below).

##
## When to consider a Bayesian approach

The York (2004) solution implemented by this package is a fast, exact generalized-least-squares fit. It is the method of choice whenever the per-point covariance matrices are trustworthy and the points scatter around the line only because of measurement error. A simulation study (correlated X-Y errors, n = 60, 40 replicates) confirmed that in this regime York and a full Bayesian treatment of the same model agree to ~10⁻³ on the slope.

A Bayesian (structural) model — a latent line with an estimated intrinsic scatter `σstruct`, fitted with MCMC — becomes important in three situations:

1. **The goodness of fit `S` is well above 1.** Two things then happen that the York solution cannot express:
   - the point estimate loses optimality: in the simulation study with intrinsic scatter (σstruct = 0.5 on a line with slope 1.5), York's slope is biased high by about +0.05 on average, while a hierarchical model that estimates the scatter is essentially unbiased (+0.01) with lower RMSE (0.061 vs 0.076);
   - the `S`-inflated `σb` (Cantrell, 2008) is a first-order correction, not a full propagation of the scatter's own uncertainty. The Bayesian posterior does this exactly: for the Amen (2012) dataset (S = 3.5), the scatter model widens the slope uncertainty from σb = 0.023 to 0.13 and shifts the slope from 1.18 to 1.06.

2. **The error correlations `r` are not well known.** Misspecifying `r` is dangerous precisely because it is silent: the goodness of fit `S` does not detect it. In a dedicated study (declared r = 0.6, actual r = 0), `S` stayed at a "normal-looking" 1.4-1.9, while the phantom correlation was silently absorbed as a fake intrinsic scatter (a constant ~0.3-0.5σ overestimate) with a systematic slope bias of about -0.15. Only a model that estimates both the scatter term and the correlation can reveal the conflict, through the visible anti-correlation of their posteriors.

3. **Downstream uncertainty propagation matters.** York reports a linearized covariance for (a, b). When the fit feeds a larger pipeline (calibration chains, hierarchical models, forecasts), posterior samples propagate the non-Gaussian tails: the σstruct direction is skewed and the (a, b) confidence corridor is curved rather than elliptical.

In short: **stay with `linearfitxy` when the covariance matrices are known and `S ≈ 1`; add a structural Bayesian layer when `S >> 1`, when `r` is uncertain, or when the fit feeds an uncertainty-propagation chain.**

Practical MCMC settings (Turing.jl ≥ 0.41): initialize the chains at the York solution (`initial_params = InitFromParams((a = st.a, b = st.b, μ = mean(X), τ = std(X)))`), use `NUTS(0.9)` for the marginalized model, and prefer the marginalized formulation of the likelihood (latent X* integrated analytically: only 4-5 parameters). With the default prior-based initialization, NUTS can wander outside the narrow posterior valley and produce diverging chains.

##
## Performance with large datasets

- The York fit itself is extremely fast: each weight-iteration pass is O(n) and converges in a few iterations regardless of n. A million points fit in well under a second; the fit is never the bottleneck.
- The Bayesian formulation matters at large n:
  - avoid the joint model that samples the n latent X*ᵢ values (O(n) parameters): beyond a few thousand points it becomes slow and mixes poorly;
  - the marginalized model has only 4-5 parameters; each likelihood evaluation is a vectorized O(n) pass, so n = 10⁵ is routine. Use the closed-form 2×2 per-point log-density instead of constructing an `MvNormal` object per point for an order-of-magnitude speedup.

##
## Note on Amen (2012)

While the results are consistent with the references provided, Amen (2012) reports regression standard deviations for Example-2 that are much smaller than this package. The discrepancy is now understood: it equals the reduced chi-square `S = 3.5` of the fit. Amen's figures correspond to the formal (uninflated) York covariance, while this package applies the Cantrell (2008) `S`-inflation by design (0.058 / 3.49 ≈ 0.017, the order of Amen's reported errors). The large `S` indicates that the stated covariances underestimate the true dispersion: an explicit intrinsic-scatter treatment of this dataset shifts the slope from 1.18 to about 1.06 (estimated σstruct ≈ 0.95).

##
## References:

*Altman, D. and Gardner, M. [1988] Statistics in Medicine: Calculating confidence intervals for regression and correlation. British Medical Journal, 296(6631), pp.1238–1242.*

*Amen, S. [2012] Linear estimation for data with error ellipses. MSc. Statistics, Univ. of Texas*

*Cantrell, C. [2008] Technical Note: Review of methods for linear least-squares fitting of data and application to atmospheric chemistry problems. Atmospheric Chem. & Physics, 8(17), pp.5477–5487*

*Mahon, K. [1996] The New "York" Regression: Application of an Improved Statistical Method to Geochemistry. International Geology Review, 38(4), pp.293–303*

*Reduced Chi-squared Test: https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic*

*Regression dilution: https://en.wikipedia.org/wiki/Regression_dilution*

*York, D. [1966] Least-squares fitting of a straight line. Canadian Journal of Physics, 44(5), pp.1079–1086*

*York, D. [1969] Least squares fitting of a straight line with correlated errors. Earth and Planetary Science Letters, 5, pp.320–324.*

*York, D., Evensen, N., Martinez, M. and Delgado J. [2004] Unified equations for the slope; intercept and standard errors of the best straight line. Am. J.Phys. 72 [3]*

##
## Example-0: no errors in X and Y
*References: pba.ucdavis.edu (0a) and Altman and Gardner (1988) (0b)*

Without input errors, the package reproduces OLS and draws the pointwise 95%
confidence band of the fitted line.

<!-- UPLOAD: Example0a_LinearFitXYerrors.png -->

![Example0b_LinearFitXYerrors](https://user-images.githubusercontent.com/20739393/132855741-7fbc7d80-76af-4f0b-a960-13e12872f8fd.png)


## Example-1: uncorrelated errors in X and Y
*References: York (1966) and Cantrell (2008)*

Data with errors on both variables and the covariance-error ellipses. The
figure also overlays the fits obtained by declaring errors in only one of the
variables (green: Y errors only; orange: X errors only), showing the
regression-dilution effect on the slope.

![Example1_LinearFitXYerrors](https://user-images.githubusercontent.com/20739393/132855764-ce6f5e40-8d5f-4d18-9c10-bdd23fec927b.png)


## Example-2: correlated errors in X and Y
*Reference: Amen (2012)*

Per-point error correlations (including negative values), covariance-error
ellipses, and the 95% confidence band of the fitted line.

Note: the shaded band is the **confidence band of the fitted line** — it
bounds the line's uncertainty, not the data's spread — *not* a prediction band
for future data points. Points may legitimately fall outside it when the stated
errors underestimate the true dispersion — which is the case here, as the large
reduced chi-square (S ≈ 3.5) signals: see
*[When to consider a Bayesian approach](#when-to-consider-a-bayesian-approach)*
for an explicit intrinsic-scatter treatment of this dataset.

![Example2_LinearFitXYerrors](https://user-images.githubusercontent.com/20739393/132855799-20aeb2fb-c327-46c6-b864-dc5ebca75736.png)


## Example-3: correlated vs uncorrelated errors
*Reference: Mahon (1996)*

Same data fitted with uncorrelated errors (r = 0):

![Example3_LinearFitXYerrors](https://user-images.githubusercontent.com/20739393/132855823-f056ce55-c017-4360-9ecf-b5fbcd6f3582.png)

and with correlated errors (r = 0.7071):

<!-- UPLOAD: Example3b_LinearFitXYerrors.png -->


## Example-4: errors in only one variable (new in v0.2.0)
*Synthetic data with known ground truth: a = 2, b = 1.5, σ = 0.4*

What happens when errors on one side are declared but not the other? Each
title displays the recovered slope versus the true one (1.5):

- errors declared in Y only (line attenuated, errors in X ignored):

<!-- UPLOAD: Example4a_LinearFitXYerrors.png -->

- errors declared in X only:

<!-- UPLOAD: Example4b_LinearFitXYerrors.png -->

- errors declared in both X and Y (the York fit recovers the true slope,
  overlaid with the fits that ignore errors in X or in Y):

<!-- UPLOAD: Example4c_LinearFitXYerrors.png -->