# LinearFitXYerrors.jl

This Julia package performs 1D linear fitting to experimental data: Y = a + b*Y\
The input data can have errors in both X and Y (error ellipses)\
The errors can be correlated and a bivariate Gaussian distribution is assumed.

This package estimates:
- The slope `b` and intercept `a` and their uncertainties
- Goodness of fit
- Pearson's correlation coefficient

The default argument `plot=false` can be turned on to plot the results.\
Currently `using Plots.jl; gr()` is used.\

## Installation
```julia
julia> ] add https://github.com/rafael-guerra-www/LinearFitXYerrors.jl
julia> using LinearFitXYerrors
```

## Useage
```julia
a, b, σa, σb, S, ρ = linearfit_xy_errors(X, Y, σX, σY)
a, b, σa, σb, S, ρ = linearfit_xy_errors(X, Y, σX, σY; r=0, plot=false)
a, b, σa, σb, S, ρ = linearfit_xy_errors(X, Y, σX, σY; r=0, plot=false)
```
#  Y = a + bX      linear fit with errors in both X & Y

#  X Y are input vectors with length > 2

#  σX and σY are vectors or scalars with the standard deviation errors in X and Y
#  Pearson's correlation coefficient r = covXY / (σX * σY)
#  r can be a vector or scalar (for constant covariance)
#  The probability distribution of each measurement is assumed to be a bivariate Gaussian

#  Ŝ is a measure of goodness of fit, perfect linear fit if Ŝ = 1
#  ρ is Pearson correlation coefficient taking in account data errors

```


## References:
*Amen, S.K. [2012] Linear estimation for data with error ellipses. MSc. Statistics, Univ. of Texas*

*Cantrell, C. [2008] Technical Note: Review of methods for linear least-squares fitting of data and application to atmospheric chemistry problems. Atmospheric Chem. & Physics, 8(17), pp.5477–5487*

*York, D. [1966] LEAST-SQUARES FITTING OF A STRAIGHT LINE. Canadian Journal of Physics, 44(5), pp.1079–1086*

*York, D., Evensen, N., Martinez, M. and Delgado J. [2004] Unified equations for the slope; intercept and standard errors of the best straight line. Am. J.Phys. 72 [3]*

*Regression dilution, https://en.wikipedia.org/wiki/Regression_dilution*

### Example-1a in examples folder produces:
![LinearFitXYerrors_example1a](https://user-images.githubusercontent.com/20739393/131935038-81db52a3-a9e5-43ab-b28b-1b701b11952f.png)
#
### Example-1b in examples folder produces:
![LinearFitXYerrors_example1b](https://user-images.githubusercontent.com/20739393/131935054-eab90824-c892-485c-9dd3-e26d61b434e7.png)
#
### Example-2 in examples folder produces:
![LinearFitXYerrors_example2](https://user-images.githubusercontent.com/20739393/131934790-68da2f2e-b132-4d65-89a6-54e92c324db2.png)
