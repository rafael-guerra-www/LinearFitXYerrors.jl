##############   EXAMPLE-0    ##############
# No errors in X nor in Y

using LinearFitXYerrors

# INPUT DATA:
# Table I from Altman and Gardner [1988]

X = [91., 104, 107, 107, 106, 100, 92, 92, 105, 108]
Y = [9.8, 7.4, 7.9, 8.3, 8.3, 9.0, 9.7, 8.8, 7.6, 6.9]

σX =  σY = r = 0


# COMPUTE and PLOT:
stxy = linearfitxy(X,Y, σX=σX, σY=σY, r=r, isplot=true)


########################################
