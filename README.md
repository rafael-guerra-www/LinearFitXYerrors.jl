# LinearFitXYerrors.jl

<img src="https://user-images.githubusercontent.com/20739393/132685450-1a34351f-ad02-49ee-9d57-b498437e356a.png" width="320" />

This small Julia package, based on York et al. (2004), performs 1D linear fitting of experimental data with uncertainties in both X and Y:

            Linear fit:             Y = a + b*X                             [1]
            
            Errors:                 X ± σX;  Y ± σY                         [2]

            Errors' correlation:    r =  = cov(σX, σY) / (σX * σY)          [3]

where:
- `X` and `Y` are input data vectors with length ≥ 3
- Optional standard deviation errors `σX` and `σY` are vectors or scalars
- Optional `r` is the correlation between the `σX` and `σY` errors.\
           `r` can be a vector or scalar

For the `σX` and `σY` errors (error ellipses) a bivariate Gaussian distribution is assumed.\
If no errors are provided, or if only `σX` or `σY` are provided, then equivalent regression results could be obtained using the [LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl) package.

The package computes:
- The intercept `a`, the slope `b` and their uncertainties `σa` and `σb`
- `σa95` and `σb95`: 95%-confidence intervals using a two-tailed t-Student distribution, e.g.: `b ± σb95 = b ± t(0.975,N-2)*σb`
- Goodness of fit `S` (reduced Χ² test): Standard Error of Estimate with Χ² distribution with N-2 degrees of freedom\
  `S ~ 1`: fit consistent with errors, `S > 1`: poor fit, `S >> 1`: errors underestimated, `S < 1`: overfitting or errors overestimated
- Pearson's correlation coefficient `ρ` that accounts for data errors
- Optional display of results with error ellipses and confidence intervals (the latter for no errors case only)

The default argument `isplot=false` can be turned on to plot the results.\
Currently `using Plots.jl; gr()` is used.


##
## Installation
```julia
julia> ] add LinearFitXYerrors
julia> using LinearFitXYerrors
```
##
## Useage
```julia
# The input data and regression results are returned in the fields of the `st` structure (::stfitxy):

st = linearfitxy(X, Y)    # no errors in X and Y, no plot displayed

st = linearfitxy(X, Y; σX, σY, isplot=true)    # X-Y errors non-correlated (r=0); plot with ratio=1

st = linearfitxy(X, Y; σX, σY, r=0, isplot=true, ratio=:auto)  # X-Y errors non-correlated (r=0); plot with auto ratio
```

## Notes:
- The objective for this first package was to learn how to publish a Julia package via Github while implementing York's technique.
- Currently the confidence interval "hyperbolic" plot ribbons are only provided for the case where input data have no errors, but in all cases, linear ribbons accounting for the standard deviation of the regression results are produced.
- The package author is not a statistician and the topics of "errors in variables" and "confidence intervals" are beyond his expertise.
- While the results seem to be consistent with the references provided, one notable exception is Amen (2012). The latter estimates standard deviation errors for regression in Example-2 that seem to be much smaller than using this package technique (York et al., 2004). However, the input data in that example have large correlated errors and York's solution seems reasonable (tbc).


##
## References:

*Altman, D. and Gardner, M. [1988] Statistics in Medicine: Calculating confidence intervals for regression and correlation. British Medical Journal, 296(6631), pp.1238–1242.*

*Amen, S. [2012] Linear estimation for data with error ellipses. MSc. Statistics, Univ. of Texas*

*Cantrell, C. [2008] Technical Note: Review of methods for linear least-squares fitting of data and application to atmospheric chemistry problems. Atmospheric Chem. & Physics, 8(17), pp.5477–5487*

*Mahon, K. [1996] The New “York” Regression: Application of an Improved Statistical Method to Geochemistry. International Geology Review, 38(4), pp.293–303*

*Reduced Chi-aquared Test: https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic*

*Regression dilution: https://en.wikipedia.org/wiki/Regression_dilution*

*York, D. [1966] Least-squares fitting of a straight line. Canadian Journal of Physics, 44(5), pp.1079–1086*

*York, D. [1969] Least squares fitting of a straight line with correlated errors. Earth and Planetary Science Letters, 5, pp.320–324.*

*York, D., Evensen, N., Martinez, M. and Delgado J. [2004] Unified equations for the slope; intercept and standard errors of the best straight line. Am. J.Phys. 72 [3]*

##
## Example-0: no errors in X and Y
*Reference: Altman and Gardner (1988)*

![Example0b_LinearFitXYerrors](https://user-images.githubusercontent.com/20739393/132855741-7fbc7d80-76af-4f0b-a960-13e12872f8fd.png)


## Example-1: uncorrelated errors in X and Y
*References: York (1966) and Cantrell (2008)*

![Example1_LinearFitXYerrors](https://user-images.githubusercontent.com/20739393/132855764-ce6f5e40-8d5f-4d18-9c10-bdd23fec927b.png)


## Example-2: correlated errors in X and Y
*Reference: Amen (2012)*

![Example2_LinearFitXYerrors](https://user-images.githubusercontent.com/20739393/132855799-20aeb2fb-c327-46c6-b864-dc5ebca75736.png)


## Example-3: correlated errors in X and Y
*Reference: Mahon (1996)*

![Example3_LinearFitXYerrors](https://user-images.githubusercontent.com/20739393/132855823-f056ce55-c017-4360-9ecf-b5fbcd6f3582.png)

