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
# The input data and regression results are returned in the fields of the `st` structure:

st = linearfitxy(X, Y)    # no errors in X and Y, no plot displayed

st = linearfitxy(X, Y; σX, σY, isplot=true)    # X-Y errors not correlated (r=0); plot with ratio=1

st = linearfitxy(X, Y; σX, σY, r=0, isplot=true, ratio=:auto)  # X-Y errors not correlated (r=0); plot with ratio=1
```

## NOTES:
- The objective for this package was to learn how to publish a Julia package via Github.
- While the results seem consistent with the references provided, one exception is Amen (2012). The latter estimates standard deviation errors that are much smaller than York et al. (2004), for input data with large errors, which are correlated. See references for further details.

- Currently confidence interval ribbons are only provided for input data with no errors.

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
## Example-0b: no errors in X and Y
*Reference: Altman and Gardner (1988)*

![Example0_LinearFitXYerrors](https://user-images.githubusercontent.com/20739393/132667125-915c4fb8-0b29-438c-a269-efeada647597.png)

## Example-1: uncorrelated errors in X and Y
*References: York (1966) and Cantrell (2008)*

![Example1_LinearFitXYerrors](https://user-images.githubusercontent.com/20739393/132667149-4cacd88a-6d62-409b-b08c-c69e78e671e3.png)

## Example-2: correlated errors in X and Y
*Reference: Amen (2012)*

![Example2_LinearFitXYerrors](https://user-images.githubusercontent.com/20739393/132667167-d4cdf29e-32a0-4e39-a990-5f00165ffc1b.png)


## Example-3: correlated errors in X and Y
*Reference: Mahon (1996)*
