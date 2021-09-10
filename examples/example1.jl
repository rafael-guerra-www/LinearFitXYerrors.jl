##############   EXAMPLE-1    ##############
# Non-correlated errors in X and in Y

using LinearFitXYerrors
using Plots; gr()

# INPUT DATA:
# York, D. [1966] Least-squares fitting of a straight line. Canadian Journal of Physics, 44(5), pp.1079–1086
# Cantrell, C. [2008] Technical Note: Review of methods for linear least-squares fitting of data and application to atmospheric chemistry problems. Atmospheric Chem. & Physics, 8(17), pp.5477–5487

# Tables I and II from York [1966]

X = [0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4]
Y = [5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5]

# the weights ωX, ωY are the inverse of the variances
σX = 1 ./ sqrt.([1000., 1000, 500, 800, 200, 80,  60, 20, 1.8, 1])
σY = 1 ./ sqrt.([1., 1.8, 4, 8, 20, 20, 70, 70, 100, 500])
r = 0*X;


# COMPUTE and PLOT:
stxy = linearfitxy(X,Y; σX=σX, σY=σY, r=r, isplot=true)



# If assuming only erros in Y or in X:
sty = linearfitxy(X,Y; σX=0, σY=σY, r=r)
stx = linearfitxy(X,Y; σX=σX, σY=0, r=r)

dX = diff([extrema(X)...])[1]/7
x1, x2 = (-dX, dX) .+ extrema(X)
xx = [x1; X; x2]
plot!(xx, sty.a .+ sty.b*xx, color=:lime, lw=0.5, label="LinearFitXY (Y errors)")
plot!(xx, stx.a .+ stx.b*xx, color=:orange, lw=0.5, label="LinearFitXY (X errors)")

@printf("LinearFitXY (Y errors): Y = (%.4f ± %.4f) + (%.4f ± %.4f)*X", sty.a, sty.σa, sty.b, sty.σb)
@printf("Pearson ρ = %.2f;  Goodness of fit = %.2f", sty.ρ, sty.S)

@printf("LinearFitXY (X errors): Y = (%.4f ± %.4f) + (%.4f ± %.4f)*X", stx.a, stx.σa, stx.b, stx.σb)
@printf("Pearson ρ = %.2f;  Goodness of fit = %.2f", stx.ρ, sty.S)


savefig("Example1_LinearFitXYerrors.png")
########################################
