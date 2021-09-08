# LinearFitXYerrors.jl

This Julia package, based on York (1966) and York et al. (2004), performs 1D linear fitting of experimental data with uncertainties in both X and Y:

            Linear fit:             Y = a + b*X                             [1]
            
            Errors:                 X ± σX;  Y ± σY                         [2]
            
            Errors' correlation:    r =  = cov(σX, σY) / (σX * σY)          [3]

where:
- `X` and `Y` are input data vectors with length ≥ 2
- Optional standard deviation errors `σX` and `σY` are vectors or scalars
- Optional `r` is the correlation between the `σX` and `σY` errors\
           `r` can be a vector or scalar

For the `σX` and `σY` errors (error ellipses) a bivariate Gaussian distribution is assumed.\
If no errors are provided, or if only `σX` or `σY` are provided, then the results are equivalent to those obtained using the [LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl) package.

`LinearFitXYerrors.jl` is based on York (1966) and York et al. (2004). See references for further details.

The package computes:
- The intercept `a`, the slope `b` and their uncertainties `σa` and `σb`
- `σa95` and `σb95`: 95%-confidence interval uncertainties corrected by two-tailed t-Student distribution, e.g.: `b ± σb95 = b ± t(0.975,N-2)*σb`
- Goodness of fit `S` (reduced Χ² test): the underlying quantity has Χ² distribution with N-2 degrees of freedom\
  `S ~ 1`: fit consistent with errors, `S > 1`: poor fit, `S >> 1`: errors underestimated, `S < 1`: errors overestimated
- Pearson's correlation coefficient `ρ` that accounts for data errors
- Optional display of fit results with error ellipses and confidence intervals

The default argument `isplot=false` can be turned on to plot the results.\
Currently `using Plots.jl; gr()` is used.

##
## Installation
```julia
julia> ] add https://github.com/rafael-guerra-www/LinearFitXYerrors.jl
julia> using LinearFitXYerrors
```
##
## Useage
```julia
# regression results are in the fields of structure st:
st = linearfitxy(X, Y)    # no errors in X and Y, no plot displayed

st = linearfitxy(X, Y; σX, σY, isplot=true)    # X-Y errors not correlateed (r=0); plot with ratio=1

st = linearfitxy(X, Y; σX, σY, r=0, isplot=true, ratio=:auto) # # X-Y errors not correlateed (r=0); plot with ratio=1
```

##
## References:

*Altman, D. and Gardner, M. [1988] Statistics in Medicine: Calculating confidence intervals for regression and correlation. British Medical Journal (Clinical research ed.), [online] 296(6631), pp.1238–1242.*

*Amen, S.K. [2012] Linear estimation for data with error ellipses. MSc. Statistics, Univ. of Texas*

*Cantrell, C. [2008] Technical Note: Review of methods for linear least-squares fitting of data and application to atmospheric chemistry problems. Atmospheric Chem. & Physics, 8(17), pp.5477–5487*

*Mahon, K. [1996] The New “York” Regression: Application of an Improved Statistical Method to Geochemistry. International Geology Review, 38(4), pp.293–303*

*Reduced Chi-aquared Test: https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic*

*Regression dilution: https://en.wikipedia.org/wiki/Regression_dilution*

*Titterington, D. and Halliday, A. [1979] On the fitting of parallel isochrons and the method of maximum likelihood. Chemical Geology, 26(3), pp.183-195*

*York, D. [1966] Least-squares fitting of a straight line. Canadian Journal of Physics, 44(5), pp.1079–1086*

*York, D., Evensen, N., Martinez, M. and Delgado J. [2004] Unified equations for the slope; intercept and standard errors of the best straight line. Am. J.Phys. 72 [3]*

##
## Example-1: uncorrelated errors
![LinearFitXYerrors_example1b](https://user-images.githubusercontent.com/20739393/131935054-eab90824-c892-485c-9dd3-e26d61b434e7.png)
##
## Example-2: correlated errors
![LinearFitXYerrors_example2](https://user-images.githubusercontent.com/20739393/131934790-68da2f2e-b132-4d65-89a6-54e92c324db2.png)
