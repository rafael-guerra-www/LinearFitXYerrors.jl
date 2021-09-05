# linearfit_xy_errors.jl


struct lfit
    a  :: Float64
    b  :: Float64
    σa :: Float64
    σb :: Float64
    S  :: Float64
    ρ  :: Float64
    bᵢ :: Vector{Float64}
    i  :: Int64
  end


# TODO: confidence_intervals

"""
# a, b, σa, σb, S, ρ, bᵢ, i = linearfit_xy_errors(X,Y,σX,σY; r=0)
#
#       Y = a + bX     : linear fit with errors in both X & Y
#
# X Y are input vectors with length > 2
#
# σX and σY are vectors or scalars with the standard deviation errors in X and Y
# Pearson's correlation coefficient r = covXY / (σX * σY)
# r can be a vector or scalar (for constant covariance)
# The probability distribution of each measurement is assumed to be a bivariate Gaussian
#
# Ŝ is a measure of goodness of fit, perfect linear fit if Ŝ = 1
# ρ is Pearson correlation coefficient taking in account data errors
# bᵢ, i : to QC convergence of b and number of iterations i required
#
#
# References:
#
# Amen, S.K. [2012] Linear estimation for data with error ellipses. MSc. Statistics, Univ. of Texas
#
# Cantrell, C. [2008] Technical Note: Review of methods for linear least-squares fitting of data and
# application to atmospheric chemistry problems. Atmospheric Chem. & Physics, 8(17), pp.5477–5487
#
# York, D. [1966] Least-squares fitting of a straight line. Canadian Journal of Physics, 44(5), pp.1079–1086
#
# York, D., Evensen, N., Martinez, M. and Delgado J. [2004] Unified equations for the slope; intercept
# and standard errors of the best straight line. Am. J.Phys. 72 [3]
#
# Regression dilution, https://en.wikipedia.org/wiki/Regression_dilution
#
#
"""
function linearfit_xy_errors(X::AbstractArray{<:Real}, Y::AbstractArray{<:Real}; σX=0, σY=0, r=0, isplot=false)

    Nmax = 50;                  # maximum number of iterations
    tol = 1e-15;                #relative tolerance to stop at
    N =   length(X)
    length(σX) == 1 && (σX = (1e-20 + σX)*ones(N))
    length(σY) == 1 && (σY = (1e-20 + σY)*ones(N))
    length(r) == 1 && (r = r*ones(N))
    
    b = Y\X                     # initial guess for b:

    bᵢ = zeros(Nmax+1);         #vector to save b iterations in
    bᵢ[1] = b
    X̄ = 0; Ȳ = 0;
    β = zeros(N)
    W = zeros(N)
    local i
    for outer i = 1:Nmax        # KEY TRICK
        ωX = 1 ./σX.^2
        ωY = 1 ./σY.^2
        α = sqrt.(ωX .* ωY)
        W .= ωX .* ωY ./(ωX + b^2 * ωY - 2*b * r .* α)
        X̄ = sum(W .* X)/sum(W)
        Ȳ = sum(W .* Y)/sum(W)
        U = X .- X̄
        V = Y .- Ȳ
        β .= W .*(U ./ ωY + b * V ./ ωX - (b * U + V) .* r ./ α)
        b = sum(W .* β .* V)/sum(W .* β .* U)
        bᵢ[i+1] = b
        abs((bᵢ[i+1] - bᵢ[i])/bᵢ[i+1]) < tol && break
    end
    a = Ȳ - b*X̄
    x = X̄ .+ β                  # y = Ȳ + b*β
    X̄ = sum(W .* x)/sum(W)      # Ȳ = sum(W .* y)/sum(W)
    u = x .- X̄                  # v = y .- Ȳ    
    Ŝ = sqrt(sum(W .* (Y - b*X .- a).^2) /(N-2))  # goodness of fit
    σb = sqrt(1/sum(W .* u.^2))
    σa = sqrt(1/sum(W) + X̄^2 * σb^2) 

    # See Wikipedia Regression dilution (Pearson's correlation coefficient with errors in variables)
    vX = var(X); vY = var(Y)
    ρ = cov(X,Y)/sqrt(vX*vY) * sqrt(vX/(vX + var(σX))) * sqrt(vY/(vY + var(σY))) 


    # standard error *Ŝ factor defined in Cantrell (2008)
    return lfit(a, b, σa*Ŝ, σb*Ŝ, Ŝ, ρ, bᵢ, i)
end


function plot_covariance_ellipses!(X,Y,a, b, c; lc=:lightgrey, lw=0.3)
    # https://cookierobotics.com/007/
    # covariance matrix = [a b; b c]
    t = LinRange(0, 2π, 72)
    st = sin.(t)
    ct = cos.(t)
    for (X, Y, a, b, c) in zip(X, Y, a, b, c)
        r = sqrt(((a - c)/2)^2 + b^2)
        λ₁ = sqrt(abs((a + c)/2 + r))
        λ₂ = sqrt(abs((a + c)/2 - r))

        (b == 0 && a >= c) ? (θ = 0.) :
        (b == 0 && a < c)  ? (θ = π/2) :
        (θ = atan(λ₁^2 - a, b))

        xₜ = @. X + λ₁*cos(θ)*ct - λ₂*sin(θ)*st
        yₜ = @. Y + λ₁*sin(θ)*ct + λ₂*cos(θ)*st
        plot!(xₜ, yₜ, lc=lc, lw=0.3, label=false)
    end
    display(plot!())
end



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

