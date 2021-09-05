# LinearFitXYerrors.jl

This Julia package, based on York (1966) and York et al. (2004), performs 1D linear fitting of experimental data with uncertainties in both X and Y:

            Linear fit:             Y = a + b*X                         [1]
            
            Errors:                 X ± σX;  Y ± σY                     [2]
            
            Errors correlation:     r =  = covXY / (σX * σY)            [3]

where:
- `X` and `Y` are input data vectors with length ≥ 2
- Optional standard deviation errors `σX` and `σY` are vectors or scalars
- Optional `r` is the correlation between the `σX` and `σY` errors\
           `r` can be a vector or scalar (for constant covariance)



The `σX` and `σY` errors (error ellipses) can be correlated, a bivariate Gaussian distribution is assumed.\
If no errors are provided, or if only `σX` or `σY` are provided, then the results are equivalent to those obtained using the [LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl) package.

`LinearFitXYerrors.jl` is based on York (1966) and York et al. (2004). See references for further details.

The package delivers:
- The intercept `a`, the slope `b` and their uncertainties
- Goodness of fit `Ŝ`, perfect linear fit if `Ŝ = 1`
- Pearson's correlation coefficient `ρ` that accounts for data errors
- Plot recipe to display fit results with error ellipses

The default argument `plot=false` can be turned on to plot the results.\
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
a, b, σa, σb, S, ρ = linearfit_xy_errors(X, Y, σX, σY)
a, b, σa, σb, S, ρ = linearfit_xy_errors(X, Y, σX, σY; r=0, plot=false)
a, b, σa, σb, S, ρ = linearfit_xy_errors(X, Y, σX, σY; r=0, plot=false)
```

##
## References:
*Amen, S.K. [2012] Linear estimation for data with error ellipses. MSc. Statistics, Univ. of Texas*

*Cantrell, C. [2008] Technical Note: Review of methods for linear least-squares fitting of data and application to atmospheric chemistry problems. Atmospheric Chem. & Physics, 8(17), pp.5477–5487*

*Regression dilution, https://en.wikipedia.org/wiki/Regression_dilution*

*York, D. [1966] LEAST-SQUARES FITTING OF A STRAIGHT LINE. Canadian Journal of Physics, 44(5), pp.1079–1086*

*York, D., Evensen, N., Martinez, M. and Delgado J. [2004] Unified equations for the slope; intercept and standard errors of the best straight line. Am. J.Phys. 72 [3]*

##
## Example-1: uncorrelated errors
![LinearFitXYerrors_example1b](https://user-images.githubusercontent.com/20739393/131935054-eab90824-c892-485c-9dd3-e26d61b434e7.png)
##
## Example-2: correlated errors
![LinearFitXYerrors_example2](https://user-images.githubusercontent.com/20739393/131934790-68da2f2e-b132-4d65-89a6-54e92c324db2.png)
