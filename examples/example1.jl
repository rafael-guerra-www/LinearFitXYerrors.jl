
##############   EXAMPLE-1a     ##############
# INPUT DATA:
# Cantrell, C. [2008] Technical Note: Review of methods for linear least-squares fitting of data and
# application to atmospheric chemistry problems. Atmospheric Chem. & Physics, 8(17), pp.5477–5487
# Data taken from Pearson (1901)

X = [0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4]
Y = [5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5]
σX = 1 ./ sqrt.([1000., 1000, 500, 800, 200, 80,  60, 20, 1.8, 1])
σY = 1 ./ sqrt.([1., 1.8, 4, 8, 20, 20, 70, 70, 100, 500])
r = 0*X;


# COMPUTE:
a, b, σa, σb, Ŝ, ρ, bᵢ, ni = linearfit_xy_errors(X,Y; σX=σX, σY=σY, r=r)


# PLOT:
using Printf, Measures, Plots; gr() 

plot_font = "Computer Modern";
default(fontfamily=plot_font,framestyle=:axes,yminorgrid=true, legendtitlefontsize=6,fg_color_legend=nothing,
    legendfontsize=6, guidefont=(7,:black),tickfont=(6,:black),size=(600,400),dpi=300, margin=0mm,
    titlefont = (8, plot_font))

# plot bᵢ with convergence of slope parameter b:
plot(1:ni+1, bᵢ[1:ni+1], xlabel="Iteration #", ylabel="b - slope")

# sort data because of ribbon:
I = sortperm(X);  X = X[I];  Y = Y[I];

str1 = @sprintf("Y = (%.4f +/- %.4f) + (%.4f +/- %.4f)*X", a,  σa, b, σb)
str2 = @sprintf("\nPearson r= %.2f; Goodness of fit = %.2f", ρ, Ŝ)

x1, x2 = (-1., 1.5) .+ extrema(X)
xx = [x1; X; x2]
tl, bl = (a - σa) .+ (b + σb)*xx,   (a + σa) .+ (b - σb)*xx

σp, σm = maximum([tl bl], dims=2) .-  (a .+ b*xx),  (a .+ b*xx) .- minimum([tl bl], dims=2)

plot(xlims= (x1,x2), title=str1*str2, legend=false)
plot!(xx, a .+ b*xx, color=:lightblue, ribbon=(σp,σm))
plot!(xx, a .+ b*xx, color=:blues, lw=1, xlabel="X", ylabel="Y")
scatter!(X, Y, ms=2,  mc=:blue, xerror=σX, yerror=σY)

@printf("LinearFitXY (X-Y errors): Y = (%.4f +/- %.4f) + (%.4f +/- %.4f)*X", a, σa, b, σb)
@printf("Pearson r = %.2f; Goodness of fit = %.2f", ρ, Ŝ)




##############   EXAMPLE-1b    ##############
# INPUT DATA: Tables I and II from York [1966]
# Data taken from Pearson (1901)

X = [0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4]
Y = [5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5]

# the weights ωX, ωY are the inverse of the variances
σX = 1 ./ sqrt.([1000., 1000, 500, 800, 200, 80,  60, 20, 1.8, 1])
σY = 1 ./ sqrt.([1., 1.8, 4, 8, 20, 20, 70, 70, 100, 500])
r = 0*X;


# COMPUTE:
a, b, σa, σb, Ŝ, ρ, bᵢ, ni = linearfit_xy_errors(X,Y; σX=σX, σY=σY, r=r)


# PLOT:
using Printf, Measures, Plots; gr() 

plot_font = "Computer Modern";
default(fontfamily=plot_font,framestyle=:axes,yminorgrid=true, legendtitlefontsize=6,fg_color_legend=nothing,
    legendfontsize=6, guidefont=(7,:black),tickfont=(6,:black),size=(600,400),dpi=300, margin=0mm,
    titlefont = (8, plot_font))

# plot bᵢ with convergence of slope parameter b:
plot(1:ni+1, bᵢ[1:ni+1], xlabel="Iteration #", ylabel="b - slope")

# sort data because of ribbon:
I = sortperm(X);  X = X[I];  Y = Y[I];

str1 = @sprintf("LinearFitXY (X-Y errors): Y = (%.4f +/- %.4f) + (%.4f +/- %.4f)*X", a, σa, b, σb)
str2 = @sprintf("\nPearson r = %.2f; Goodness of fit = %.2f", ρ, Ŝ)
x1, x2 = (-1., 1.5) .+ extrema(X)
xx = [x1; X; x2]
tl, bl = (a - σa) .+ (b + σb)*xx,   (a + σa) .+ (b - σb)*xx

σp, σm = maximum([tl bl], dims=2) .-  (a .+ b*xx),  (a .+ b*xx) .- minimum([tl bl], dims=2)

plot(xlims=(x1,x2), title=str1*str2, legend=:bottomleft)
plot!(xx, a .+ b*xx, color=:lightblue, ribbon=(σp,σm), label=false)
plot!(xx, a .+ b*xx, color=:blue, lw=1, xlabel="X", ylabel="Y", label="LinearFitXY (X-Y errors)")
scatter!(X, Y, ms=2,  mc=:blue, xerror=σX, yerror=σY, label=false)

@printf("LinearFitXY (X-Y errors): Y = (%.4f +/- %.4f) + (%.4f +/- %.4f)*X", a,  σa, b, σb)
@printf("Pearson r = %.2f; Goodness of fit = %.2f", ρ, Ŝ)


# If assuming only erros in Y or in X:
ay, by, σay, σby, Ŝy, ρy = linearfit_xy_errors(X,Y; σX=0, σY=σY, r=r)
ax, bx, σax, σbx, Ŝx, ρx = linearfit_xy_errors(X,Y; σX=σX, σY=0, r=r)

plot!(xx, ay .+ by*xx, color=:lime, lw=0.5, label="LinearFitXY (Y errors)")
plot!(xx, ax .+ bx*xx, color=:orange, lw=0.5, label="LinearFitXY (X errors)")

@printf("LinearFitXY (Y errors): Y = (%.4f +/- %.4f) + (%.4f +/- %.4f)*X", ay,  σay, by, σby)
@printf("Pearson r = %.2f; Goodness of fit = %.2f", ρy, Ŝy)

@printf("LinearFitXY (X errors): Y = (%.4f +/- %.4f) + (%.4f +/- %.4f)*X", ax,  σax, bx, σbx)
@printf("Pearson r = %.2f; Goodness of fit = %.2f", ρx, Ŝx)


########################################