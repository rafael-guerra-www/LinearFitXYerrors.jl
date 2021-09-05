
##############   EXAMPLE-2     ##############

using LinearFitXYerrors

# INPUT DATA: Table 3.1 from Amen, S. [2012]
# download csv file from examples folder to your local disk and edit path:

using CSV, DataFrames
file = raw".\examples\Amen_Errors_X_Y_and_covariance_Table_31.csv"
M = CSV.read(file, DataFrame)

μₑ = M[:,:ue]
η = M[:,:n]
σμₑ = sqrt.(M[:,:var_ue])   # standard deviation of errors in μₑ
ση = sqrt.(M[:,:var_n])     # standard deviation of errors in η
cov_μₑ_η = M[:,:cov_ue_n]   # covariance between errors in μₑ and η
r = M[:,:cov_ue_n] ./ (σμₑ .* ση )   # correlation coefficient between errors


# COMPUTE:
stxy = linearfit_xy_errors(μₑ, η; σX=σμₑ, σY=ση, r=r, isplot=true)


# PLOT:
plotlinfitxy(μₑ, η; σX=σμₑ, σY=ση, r=r, st=stxy)

# If assuming only erros in Y or in X:
sty = linearfit_xy_errors(μₑ, η; σX=0, σY=ση, r=r)
stx = linearfit_xy_errors(μₑ, η; σX=σμₑ, σY=0, r=r)

dX = diff([extrema(μₑ)...])[1]/7
x1, x2 = (-dX, dX) .+ extrema(μₑ)
xx = [x1; X; x2]
plot!(xx, sty.a .+ sty.b*xx, color=:lime, lw=0.5, label="LinearFitXY (Y errors)");
plot!(xx, stx.a .+ stx.b*xx, color=:orange, lw=0.5, label="LinearFitXY (X errors)");

@printf("LinearFitXY (Y errors): Y = (%.4f +/- %.4f) + (%.4f +/- %.4f)*X", sty.a, sty.σa, sty.b, sty.σb);
@printf("Pearson r = %.2f; Goodness of fit = %.2f", sty.ρ, sty.S)

@printf("LinearFitXY (X errors): Y = (%.4f +/- %.4f) + (%.4f +/- %.4f)*X", stx.a, stx.σa, stx.b, stx.σb);
@printf("Pearson r = %.2f; Goodness of fit = %.2f", stx.ρ, stx.S)


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
