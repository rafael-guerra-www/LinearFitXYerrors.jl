##############   EXAMPLE-0    ##############
# No errors in X nor in Y

using LinearFitXYerrors

# INPUT DATA:
# Altman, D. and Gardner, M. [1988] Statistics in Medicine: Calculating confidence intervals for regression and correlation. British Medical Journal, 296(6631), pp.1238–1242.

# Table I
X = [91., 104, 107, 107, 106, 100, 92, 92, 105, 108]
Y = [9.8, 7.4, 7.9, 8.3, 8.3, 9.0, 9.7, 8.8, 7.6, 6.9]


# COMPUTE and PLOT:
stxy = linearfitxy(X,Y, isplot=true)    # σX =  σY = r = 0

# savefig("Example0_LinearFitXYerrors.png")
########################################
