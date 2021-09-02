module LinearFitXYerrors

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
# York, D. [1966] LEAST-SQUARES FITTING OF A STRAIGHT LINE. Canadian Journal of Physics, 44(5), pp.1079–1086
#
# York, D., Evensen, N., Martinez, M. and Delgado J. [2004] Unified equations for the slope; intercept
# and standard errors of the best straight line. Am. J.Phys. 72 [3]
#
# Regression dilution, https://en.wikipedia.org/wiki/Regression_dilution
#
#
"""
function linearfit_xy_errors(X,Y; σX=0, σY=0, r=0)
    Nmax = 50; # maximum number of iterations
    tol = 1e-15;    #relative tolerance to stop at
    N =   length(X)
    length(σX) == 1 && (σX = (1e-20 + σX)*ones(N))
    length(σY) == 1 && (σY = (1e-20 + σY)*ones(N))
    length(r) == 1 && (r = r*ones(N))
    
    b = Y\X    # initial guess for b:

    bᵢ = zeros(Nmax+1);    #vector to save b iterations in
    bᵢ[1] = b
    X̄ = 0; Ȳ = 0;
    β = zeros(N)
    W = zeros(N)
    local i
    for outer i = 1:Nmax   # KEY TRICK
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

    # standard error *Ŝ factor defined in Cantrell
    return a, b, σa*Ŝ, σb*Ŝ, Ŝ, ρ, bᵢ, i
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

end
