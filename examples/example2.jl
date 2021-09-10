##############   EXAMPLE-2     ##############
# Correlated errors in X and Y

using LinearFitXYerrors
using Plots; gr()


# INPUT DATA:
# Amen, S. [2012] Linear estimation for data with error ellipses. MSc. Statistics, Univ. of Texas

#  Table 3.1

# Num,    ue,    n,   var_ue,   var_n,   cov_ue_n
M = [1   23.371   5.992   0.106   0.323   0.170;
    2    24.462   8.899   1.692   10.195  3.887;
    3    22.414   5.361   0.124   0.460   0.207;
    4    27.628   10.271  1.088   1.460   1.255;
    5    23.802   7.114   0.107   0.496   0.191;
    6    23.092   7.975   0.281   1.651   0.582;
    7    23.263   7.486   0.316   1.695   0.661;
    8    20.843   4.169   0.389   0.335   0.357;
    9    24.421   5.645   0.381   0.716   0.507;
    10   23.716   9.220   0.733   3.656   1.401;
    11   22.032   5.358   0.094   0.268   0.106;
    12   21.448   3.166   0.167   0.170   0.001;
    13   21.631   4.002   0.127   0.106   0.082;
    14   19.947   2.073   0.015   0.015   0.001;
    15   20.652   3.338   0.040   0.104   0.021;
    16   20.588   2.028   0.007   0.011   0.002;
    17   20.745   1.978   0.010   0.010  -0.001;
    18   21.658   2.526   0.009   0.058  -0.003;
    19   19.539   2.043   0.044   0.070  -0.039;
    20   21.638   3.361   1.909   3.961  -2.557;
    21   19.863   2.452   0.014   0.019  -0.010;
    22   18.398   2.199   0.042   0.051  -0.037;
    23   20.729   1.725   0.013   0.026   0.003;
    24   20.511   1.906   0.003   0.010   0.001;
    25   22.140   3.372   0.047   0.108   0.062;
    26   20.510   2.130   0.004   0.020   0.002;
    27   20.141   1.897   0.006   0.011  -0.005]

μₑ = M[:,2]
η = M[:,3]
σμₑ = sqrt.(M[:,4])   # standard deviation of errors in μₑ
ση = sqrt.(M[:,5])     # standard deviation of errors in η
cov_μₑ_η = M[:,6]   # covariance between errors in μₑ and η
r = cov_μₑ_η ./ (σμₑ .* ση )   # correlation coefficient between errors


# COMPUTE and PLOT:
# The standard deviation errors using York (2004) are >> than Amen (2012)
# However, Amen (2012) errors seem to be underestimated - tbc

stxy = linearfitxy(μₑ, η; σX=σμₑ, σY=ση, r=r, isplot=true)



# If assuming no errors:
st = linearfitxy(μₑ, η)

# If assuming only errors in Y or in X:
sty = linearfitxy(μₑ, η; σX=0, σY=ση, r=r)
stx = linearfitxy(μₑ, η; σX=σμₑ, σY=0, r=r)

dX = diff([extrema(μₑ)...])[1]/7
x1, x2 = (-dX, dX) .+ extrema(μₑ)
xx = [x1; μₑ; x2]
plot!(xx, st.a .+ st.b*xx, color=:pink, lw=0.5, label="LinearFitXY (no errors)");
plot!(xx, sty.a .+ sty.b*xx, color=:lime, lw=0.5, label="LinearFitXY (Y errors)");
plot!(xx, stx.a .+ stx.b*xx, color=:orange, lw=0.5, label="LinearFitXY (X errors)");

@printf("LinearFitXY (Y errors): Y = (%.4f ± %.4f) + (%.4f ± %.4f)*X", sty.a, sty.σa, sty.b, sty.σb);
@printf("Pearson r = %.2f; Goodness of fit = %.2f", sty.ρ, sty.S)

@printf("LinearFitXY (X errors): Y = (%.4f ± %.4f) + (%.4f ± %.4f)*X", stx.a, stx.σa, stx.b, stx.σb);
@printf("Pearson r = %.2f; Goodness of fit = %.2f", stx.ρ, stx.S)


# Compare with LsqFit:
using LsqFit
m(x, p) = p[1] .+ p[2] * x
p0 = [0., 0.]
fit2 = curve_fit(m, μₑ, η, p0)
cf2 = coef(fit2)
ci2 = confidence_interval(fit2, 0.05)    # 5% significance level
@printf("LsqFit (no errors): Y = (%.2f ± %.2f) + (%.2f ± %.2f)*X",
          cf2[1],diff([ci2[1]...])[1]/2, cf2[2],diff([ci2[2]...])[1]/2)

wt3 = 1 ./ ση .^2
fit3 = curve_fit(m, μₑ, η, wt3, p0)
cf3 = coef(fit3)
ci3 = confidence_interval(fit3, 0.05)    # 5% significance level
@printf("LsqFit (Y-errors): Y = (%.2f ± %.2f) + (%.2f ± %.2f)*X",
          cf3[1],diff([ci3[1]...])[1]/2, cf3[2],diff([ci3[2]...])[1]/2)

wt4 = 1 ./ σμₑ .^2
fit4 = curve_fit(m, η, μₑ, wt4, p0)     # fit X function of Y
cf4 = coef(fit4)
cf4[1] = -cf4[1] / cf4[2]       # convert coefficients to Y = a + b*X form
cf4[2] = 1 / cf4[2] 
ci4 = confidence_interval(fit4, 0.05)    # 5% significance level
@printf("LsqFit (X-errors): Y = (%.2f ± %.2f) + (%.2f ± %.2f)*X",
          cf4[1],diff([ci4[1]...])[1]/2, cf4[2],diff([ci4[2]...])[1]/2)

plot!(xx, cf2[1] .+ cf2[2]*xx, lw=0.5, lc=:red, ls=:dash, label="LsqFit (no errors)")
plot!(xx, cf3[1] .+ cf3[2]*xx, color=:darkgreen, lw=0.5, ls=:dash, label="LsqFit (Y-errors)")
plot!(xx, cf4[1] .+ cf4[2]*xx, color=:darkred, lw=0.5, ls=:dash, label="LsqFit (X-errors)")


savefig("Example2_LinearFitXYerrors.png")
########################################
