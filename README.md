# LinearFitXYerrors.jl
#
### Example-1 in examples folder produces:
![LinearFitXYerrors_example1](https://user-images.githubusercontent.com/20739393/131933542-9927aacb-32e1-433e-8fdf-2002896dfdf0.png)
#
### Example-2 in examples folder produces:
![LinearFitXYerrors_example2](https://user-images.githubusercontent.com/20739393/131933586-4739bf00-7cd8-4aaf-8e40-24626b989d84.png)

##### a, b, σa, σb, S, ρ, bᵢ, i = linearfit_xy_errors(X,Y,σX,σY; r=0)
#####
#####       Y = a + bX     : linear fit with errors in both X & Y
#####
##### X Y are input vectors with length > 2
#####
##### σX and σY are vectors or scalars with the standard deviation errors in X and Y
##### Pearson's correlation coefficient r = covXY / (σX * σY)
##### r can be a vector or scalar (for constant covariance)
##### The probability distribution of each measurement is assumed to be a bivariate Gaussian
#####
##### Ŝ is a measure of goodness of fit, perfect linear fit if Ŝ = 1
##### ρ is Pearson correlation coefficient taking in account data errors
##### bᵢ, i : to QC convergence of b and number of iterations i required
#####
#####
### References:
#####
##### Amen, S.K. [2012] Linear estimation for data with error ellipses. MSc. Statistics, Univ. of Texas
#####
##### Cantrell, C. [2008] Technical Note: Review of methods for linear least-squares fitting of data and
##### application to atmospheric chemistry problems. Atmospheric Chem. & Physics, 8(17), pp.5477–5487
#####
##### York, D. [1966] LEAST-SQUARES FITTING OF A STRAIGHT LINE. Canadian Journal of Physics, 44(5), pp.1079–1086
#####
##### York, D., Evensen, N., Martinez, M. and Delgado J. [2004] Unified equations for the slope; intercept
#### and standard errors of the best straight line. Am. J.Phys. 72 [3]
#####
##### Regression dilution, https://en.wikipedia.org/wiki/Regression_dilution
#
