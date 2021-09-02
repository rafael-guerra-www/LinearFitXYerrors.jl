
##############   EXAMPLE-2     ##############
# INPUT DATA: Table 3.1 from Amen, S. [2012]

using CSV, DataFrames

file = raw"C:\Users\jrafa\OneDrive\Julia_Code\TO_DO\Linear_fit__errors_x_y\Amen_Errors_X_Y_and_covariance_Table_31.csv"
M = CSV.read(file, DataFrame)

μₑ = M[:,:ue]
η = M[:,:n]
σμₑ = sqrt.(M[:,:var_ue])
ση = sqrt.(M[:,:var_n])
cov_μₑ_η = M[:,:cov_ue_n]
r = M[:,:cov_ue_n] ./ (σμₑ .* ση )   # Pearson's correlation coefficient


# COMPUTE:
a, b, σa, σb, Ŝ, ρ, bᵢ, ni = linearfit_xy_errors(μₑ, η; σX=σμₑ, σY=ση, r=r)


# PLOT:
using Printf, Measures, Plots; gr(dpi=300) 

plot_font = "Computer Modern";
default(fontfamily=plot_font,framestyle=:axes,yminorgrid=true, legendtitlefontsize=6,fg_color_legend=nothing,
    legendfontsize=6, guidefont=(7,:black),tickfont=(6,:black),size=(600,400),dpi=300, margin=0mm,
    titlefont = (6, plot_font))

# plot bᵢ with convergence of slope parameter b:
plot(1:ni+1, bᵢ[1:ni+1], xlabel="Iteration #", ylabel="b - slope")

# sort data because of ribbon:
I = sortperm(μₑ);
μₑ = μₑ[I];  η = η[I]; σμₑ = σμₑ[I];  ση = ση[I]; cov_μₑ_η = cov_μₑ_η[I];

str1 = @sprintf("LinearFitXY (X-Y errors): η = (%.4f +/- %.4f) + (%.4f +/- %.4f)*X", a, σa, b, σb)
str2 = @sprintf("\nPearson r = %.2f; Goodness of fit = %.2f", ρ, Ŝ)

x1, x2 = (-1., 1.5) .+ extrema(μₑ)
xx = [x1; μₑ; x2]
tl, bl = (a - σa) .+ (b + σb)*xx,   (a + σa) .+ (b - σb)*xx
σp, σm = maximum([tl bl], dims=2) .-  (a .+ b*xx),  (a .+ b*xx) .- minimum([tl bl], dims=2)


plot(xlims=(x1,x2), ylims=(0,13), title=str1*str2, ratio=1, legend=:outerbottomright)
plot!(xx, a .+ b*xx, color=:lightblue, ribbon=(σp,σm), label=false)
plot!(xx, a .+ b*xx, color=:blue, lw=0.5, xlabel="μe", ylabel="η", label="LinearFitXY")
scatter!(μₑ, η, msw=0.1, ms=1., msc=:lightgrey, xerror= σμₑ, yerror= ση, label=false)
scatter!(μₑ, η, msw=0.1, ms=1.5, mc=:blue, label=false)
plot_covariance_ellipses!(μₑ, η, σμₑ.^2, cov_μₑ_η, ση.^2; lc=:grey)

@printf("LinearFitXY (X-Y errors): η = (%.4f +/- %.4f) + (%.4f +/- %.4f)*μe", a,  σa, b, σb)
@printf("Pearson r = %.2f; Goodness of fit = %.2f", ρ, Ŝ)

# If assuming only erros in Y or in X:
ay, by, σay, σby, Ŝy, ρy = linearfit_xy_errors(μₑ, η; σX=0, σY=ση, r=r);
ax, bx, σax, σbx, Ŝx, ρx = linearfit_xy_errors(μₑ, η; σX=σμₑ, σY=0, r=r)

plot!(xx, ay .+ by*xx, color=:lime, lw=0.5, label="LinearFitXY (Y errors)");
plot!(xx, ax .+ bx*xx, color=:orange, lw=0.5, label="LinearFitXY (X errors)");

@printf("LinearFitXY (Y errors): Y = (%.4f +/- %.4f) + (%.4f +/- %.4f)*X", ay,  σay, by, σby);
@printf("Pearson r = %.2f; Goodness of fit = %.2f", ρy, Ŝy)

@printf("LinearFitXY (X errors): Y = (%.4f +/- %.4f) + (%.4f +/- %.4f)*X", ax,  σax, bx, σbx);
@printf("Pearson r = %.2f; Goodness of fit = %.2f", ρx, Ŝx)


# Compare with LsqFit:
using LsqFit
m(x, p) = p[1] .+ p[2] * x
p0 = [0., 0.]
fit2 = curve_fit(m, μₑ, η, p0)
cf2 = coef(fit2)
ci2 = confidence_interval(fit2, 0.05)    # 5% significance level
@printf("LsqFit (no errors): Y = (%.2f +/- %.2f) + (%.2f +/- %.2f)*X",
          cf2[1],diff([ci2[1]...])[1]/2, cf2[2],diff([ci2[2]...])[1]/2)

wt3 = 1 ./ ση .^2
fit3 = curve_fit(m, μₑ, η, wt3, p0)
cf3 = coef(fit3)
ci3 = confidence_interval(fit3, 0.05)    # 5% significance level
@printf("LsqFit (Y-errors): Y = (%.2f +/- %.2f) + (%.2f +/- %.2f)*X",
          cf3[1],diff([ci3[1]...])[1]/2, cf3[2],diff([ci3[2]...])[1]/2)

wt4 = 1 ./ σμₑ .^2
fit4 = curve_fit(m, η, μₑ, wt4, p0)     # fit X function of Y
cf4 = coef(fit4)
cf4[1] = -cf4[1] / cf4[2]       # convert coefficients to Y = a + b*X form
cf4[2] = 1 / cf4[2] 
ci4 = confidence_interval(fit4, 0.05)    # 5% significance level
@printf("LsqFit (X-errors): Y = (%.2f +/- %.2f) + (%.2f +/- %.2f)*X",
          cf4[1],diff([ci4[1]...])[1]/2, cf4[2],diff([ci4[2]...])[1]/2)

plot!(xx, cf2[1] .+ cf2[2]*xx, lw=0.5, lc=:red, ls=:dash, label="LsqFit (no errors)")
plot!(xx, cf3[1] .+ cf3[2]*xx, color=:darkgreen, lw=0.5, ls=:dash, label="LsqFit (Y-errors)")
plot!(xx, cf4[1] .+ cf4[2]*xx, color=:darkred, lw=0.5, ls=:dash, label="LsqFit (X-errors)")


########################################