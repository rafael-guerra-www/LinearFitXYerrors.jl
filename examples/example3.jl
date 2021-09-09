##############   EXAMPLE-3     ##############
# Correlated errors in X and Y

using LinearFitXYerrors

# INPUT DATA:
# Mahon, K. [1996] The New “York” Regression: Application of an Improved Statistical Method to Geochemistry. International Geology Review, 38(4), pp.293–303

#  Table 3

# Num    X        Y
M = [0.037   0.0080
     0.035   0.0084
     0.032   0.0100
     0.040   0.0085
     0.013   0.0270
     0.038   0.0071
     0.042   0.0043
     0.030   0.0160]

X = M[:,1]
Y = M[:,2]
σX = 0.03*X     # standard deviation of errors in X
σY = 0.10*Y     # standard deviation of errors in Y
r1 = 0  # correlation coefficient between errors
r2 = 0.7071

# COMPUTE and PLOT:
stxy1 = linearfitxy(X, Y; σX=σX, σY=σY, r=r1, isplot=true, ratio=:auto)

stxy2 = linearfitxy(X, Y; σX=σX, σY=σY, r=r2, isplot=true, ratio=:auto)


# If assuming no errors:
st = linearfitxy(X, Y)


# savefig("Example3_LinearFitXYerrors.png")
########################################
