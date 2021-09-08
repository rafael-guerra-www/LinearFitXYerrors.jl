# linearfitxy.jl


struct stfit
    a   :: Float64
    b   :: Float64
    σa  :: Float64
    σb  :: Float64
    σa95:: Float64
    σb95:: Float64
    S   :: Float64
    ρ   :: Float64
    bᵢ  :: Vector{Float64}
    i   :: Int64
  end


"""
# stfit = linearfitxy(X,Y,σX,σY; r=0)
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
# stfit.a, stfit.b
# stfit.σa, stfit.σb
# stfit.S is a measure of goodness of fit, perfect linear fit if Ŝ = 1
# stfit.ρ is Pearson correlation coefficient taking in account data errors
# stfit.bᵢ, stfit.i : to QC convergence of b and number of iterations i required
#
#
# References:
#
# Amen, S.K. [2012] Linear estimation for data with error ellipses. MSc. Statistics, Univ. of Texas
#
# Cantrell, C. [2008] Technical Note: Review of methods for linear least-squares fitting of data and
# application to atmospheric chemistry problems. Atmospheric Chem. & Physics, 8(17), pp.5477–5487
#
# Regression dilution, https://en.wikipedia.org/wiki/Regression_dilution
#
# York, D. [1966] Least-squares fitting of a straight line. Canadian Journal of Physics, 44(5), pp.1079–1086
#
# York, D., Evensen, N., Martinez, M. and Delgado J. [2004] Unified equations for the slope; intercept
# and standard errors of the best straight line. Am. J.Phys. 72 [3]
#
#
"""
function linearfitxy(X, Y; σX=0, σY=0, r=0, isplot=false, ratio=1)
    
    N =   length(X)
    d = TDist(N-2)     # t-Student distribution with N-2 degrees of freedom
    cf = quantile(d, 0.975)  # correction factor for 95% confidence intervals (two-tailed distribution)

    if (σX == 0) && (σY == 0)
        σX = σY = zeros(N)
        length(r) == 1  &&  (r = r*ones(N))
        
        M = [ones(N) X]
        (a, b) = M \ Y
        X̄ = mean(X)    # Ȳ = mean(Y)

        # standard error * Ŝ factor defined in Cantrell (2008)
        Ŝ = sqrt(sum((Y - b*X .- a).^2)/(N-2))  # goodness of fit
        σb = sqrt(1/sum((X .- X̄ ).^2))
        σa = Ŝ * sqrt(1/N + X̄^2 * σb^2)
        σb *= Ŝ 

        # Pearson's correlation coefficient:
        ρ = cov(X,Y)/sqrt(var(X) * var(Y))

        st = stfit(a, b, σa, σb, cf*σa, cf*σb, Ŝ, ρ, [b], 1)

    else
        Nmax = 50;                  # maximum number of iterations
        tol = 1e-15;                #relative tolerance to stop at
        
        σX == 0 && (σX = 1e-16)
        σY == 0 && (σY = 1e-16)
        length(σX) == 1 &&  (σX = σX*ones(N))
        length(σY) == 1 &&  (σY = σY*ones(N))
        length(r) == 1  &&  (r = r*ones(N))
        
        M = [ones(N) X]
        (a, b) = M \ Y              # initial guess for b
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

        # standard error *Ŝ factor defined in Cantrell (2008)
        Ŝ = sqrt(sum(W .* (Y - b*X .- a).^2) /(N-2))  # goodness of fit
        σb = sqrt(1/sum(W .* (x .- X̄).^2))
        σa = Ŝ * sqrt(1/sum(W) + X̄^2 * σb^2)
        σb *= Ŝ

        # See Wikipedia Regression dilution (Pearson's correlation coefficient with errors in variables)
        vX = var(X); vY = var(Y)
        ρ = cov(X,Y)/sqrt(vX*vY) * sqrt(vX/(vX + var(σX))) * sqrt(vY/(vY + var(σY))) 

        st = stfit(a, b, σa, σb, cf*σa, cf*σb, Ŝ, ρ, bᵢ, i)
    end

    if isplot
        plot_linfitxy(X, Y; σX, σY, r, st, ratio=ratio)
    end

    @printf("\n>>> ± σ: Y = (%.4f ± %.4f) + (%.4f ± %.4f)*X \n", a,  σa, b, σb)
    @printf(">>> ± 95%% CI): Y = (%.4f ± %.4f) + (%.4f ± %.4f)*X \n", a, cf*σa, b, cf*σb)
    @printf(">>> Pearson ρ = %.3f;  Goodness of fit = %.3f \n\n", ρ, Ŝ)

    return st
end


function plot_covariance_ellipses!(X, Y, a,  b, c; lc=:lightgrey, lw=0.3)
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


function  plot_linfitxy(X, Y; σX=0, σY=0, r=0, st::stfit, ratio=1)
    covXY = σX .* σY .* r
    # sort data because of ribbon:
    I = sortperm(X);
    X = X[I];  Y = Y[I]; σX = σX[I];  σY = σY[I]; covXY = covXY[I];

    str1 = @sprintf("[± σ]  Y = (%.4f ± %.4f) + (%.4f ± %.4f)*X", st.a, st.σa, st.b, st.σb)
    str2 = @sprintf("\n[± 95%% CI]  Y = (%.4f ± %.4f) + (%.4f ± %.4f)*X", st.a, st.σa95, st.b, st.σb95)
    str3 = @sprintf("\nPearson ρ = %.2f; Goodness of fit = %.2f", st.ρ, st.S)

    dX = diff([extrema(X)...])[1]/7    # extends plot x-axis by 1/7 each side
    x1, x2 = (-dX, dX) .+ extrema(X)
    xx = [x1; X; x2]
    tl, bl = (st.a - st.σa) .+ (st.b + st.σb)*xx,  (st.a + st.σa) .+ (st.b - st.σb)*xx
    σp, σm = maximum([tl bl], dims=2) .-  (st.a .+ st.b*xx),  (st.a .+ st.b*xx) .- minimum([tl bl], dims=2)

    #TODO: confirm 95%-confidence interval computation for full X-Y errors
    mx = mean(X)
    dx2 = (X .- mx).^2
    cf = st.σb95 /st.σb   # to reuse already computed t-Student quantile
    σx95 = cf * st.S * sqrt.(1/length(X) .+ dx2/sum(dx2))

    plot(xlims=(x1,x2), ylims=(0,13), title=str1*str2*str3, ratio=ratio, legend=:outerbottomright)
    scatter!([x1-dX x1-dX],[0 0], mc=[:lightblue :reds], label=["95% confidence interval" "Y=(a ± σa)+(b ± σb)*X"], marker=:rect)
    plot!([x1-dX],[0], lc=:white, label=" ")   # adds space in legend
    plot!(X, st.a .+ st.b*X, color=:lightblue, ribbon=(σx95,σx95), label=false)
    plot!(xx, st.a .+ st.b*xx, color=:reds, ribbon=(σp,σm), label=false)
    plot!(xx, st.a .+ st.b*xx, color=:blue, lw=0.5, xlabel="X", ylabel="Y", label="LinearFitXY")
    scatter!(X, Y, msw=0.1, ms=1., msc=:lightgrey, xerror= σX, yerror= σY, label=false)
    scatter!(X, Y, msw=0.1, ms=1.5, mc=:blue, label=false)
    
    if σX != 0 && σY != 0
        plot_covariance_ellipses!(X, Y, σX.^2, covXY, σY.^2; lc=:grey)
    end
end
