# linearfitxy.jl


struct stfitxy
    X   :: Vector{Float64}
    Y   :: Vector{Float64}
    σX  :: Union{Float64, Vector{Float64}}
    σY  :: Union{Float64, Vector{Float64}}
    r   :: Union{Float64, Vector{Float64}}
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
linearfitxy(X, Y; σX=0, σY=0, r=0, isplot=false, ratio=1)

Performs 1D linear fitting of experimental data with uncertainties in  X and Y:
- Linear fit:             `Y = a + b*X`                               [1]
- Errors:                 ``X ± σX;  Y ± σY``                         [2]
- Errors' correlation:    ``r =  = cov(σX, σY) / (σX * σY)``          [3]

# Arguments:
- `X` and `Y` are input data vectors with length ≥ 3
- Optional standard deviation errors ``σX`` and ``σY`` are vectors or scalars
- Optional `r` is the correlation between the ``σX`` and ``σY`` errors
`r` can be a vector or scalar

``σX`` and ``σY`` errors (error ellipses) with bivariate Gaussian distribution assumed.
If no errors, or if only ``σX`` or ``σY`` are provided, then the results are equivalent
to those from the [LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl) package.

Based on York et al. (2004) with extensions (confidence intervals, diluted corr. coeff.).

# Examples:
```julia-repl
st = linearfitxy(X, Y)    # no errors in X and Y, no plot displayed

st = linearfitxy(X, Y; σX, σY, isplot=true) # XY errors not correlated (r=0); plot ratio=1

linearfitxy(X, Y; σX,σY,r=0,isplot=true,ratio=:auto) # XY errors not correlated (r=0); plot
```

The results are in the fields of the returned st::stfitxy structure:
- The intercept `a`, the slope `b` and their uncertainties `σa` and `σb`
- ``σa95`` and ``σb95``: 95%-confidence interval using two-tailed t-Student distribution,
    e.g.: ``b ± σb95 = b ± t(0.975,N-2)*σb``
- Goodness of fit `S` (reduced ``Χ²`` test): quantity with ``Χ²`` N-2 degrees of freedom
  `S ~ 1`: fit consistent with errors, `S > 1`: poor fit, `S >> 1`: errors underestimated,
  `S < 1`: overfitting or errors overestimated
- Pearson's correlation coefficient ``ρ`` that accounts for data errors
- Optional display of fit results with error ellipses and confidence intervals

The default argument `isplot=false` can be turned on to plot the results.
Currently Plots.jl's gr() is used

"""
function linearfitxy(X, Y; σX=0, σY=0, r=0, isplot=false, ratio=1)
    
    N =   length(X)
    d = TDist(N-2)     # t-Student distribution with N-2 degrees of freedom
    cf = quantile(d, 0.975)  # correction factor for 95% confidence intervals (two-tailed distribution)

    if (σX == 0) && (σY == 0)
        σX = σY = r = zeros(N)
        
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

        st = stfitxy(X, Y, σX, σY, r, a, b, σa, σb, cf*σa, cf*σb, Ŝ, ρ, [b], 1)

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

        
        # compare to: Ŝ = sqrt(sum((Y - b*X .- a).^2)/(N-2))
        Ŝ = sqrt(sum(W .* (Y - b*X .- a).^2) /(N-2))  # goodness of fit (York, 2004), (Cantrell, 2008)
               
        σb = sqrt(1/sum(W .* (x .- X̄).^2))
        σa = Ŝ * sqrt(1/sum(W) + X̄^2 * σb^2)
        σb *= Ŝ             # standard error * Ŝ factor (Cantrell, 2008)
        
        # See Wikipedia Regression dilution (Pearson's correlation coefficient with errors in variables)
        vX = var(X); vY = var(Y)
        ρ = cov(X,Y)/sqrt(vX*vY) * sqrt(vX/(vX + var(σX))) * sqrt(vY/(vY + var(σY))) 

        st = stfitxy(X, Y, σX, σY, r, a, b, σa, σb, cf*σa, cf*σb, Ŝ, ρ, bᵢ, i)
    end

    if isplot
        plot_linfitxy(st; ratio=ratio)
    end

    @printf("\n>>> [± σ]  Y = (%.4f ± %.4f) + (%.4f ± %.4f)*X \n", a,  σa, b, σb)
    @printf(">>> [± 95%% CI]  Y = (%.4f ± %.4f) + (%.4f ± %.4f)*X \n", a, cf*σa, b, cf*σb)
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


function  plot_linfitxy(st::stfitxy; ratio=1)
# set plotting defaults
    plot_font = "Computer Modern";
    default(fontfamily=plot_font,framestyle=:axes,yminorgrid=true, legendtitlefontsize=6,fg_color_legend=nothing,
            legendfontsize=6, guidefont=(7,:black),tickfont=(6,:black),size=(600,400), dpi=300,
            titlefont = (6, plot_font))
    
    X = st.X;  Y = st.Y;  σX = st.σX;  σY = st.σY;  r = st.r;
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

    #TODO: 95%-confidence interval computation for X-Y errors
    N = length(X)
    d = TDist(N-2)     # t-Student distribution with N-2 degrees of freedom
    cf = quantile(d, 0.975)  # correction factor for 95% confidence intervals (two-tailed distribution)
    σx95 = cf * st.S * sqrt.(1/N .+ (X .- mean(X)).^2 / var(X) /(N-1))
    
    
    plot(xlims=(x1,x2), title=str1*str2*str3, ratio=ratio, legend=:outerbottomright)
    scatter!([x1-dX x1-dX],[Y[1] Y[1]], mc=[:lightblue :reds], label=["95% confidence interval" "Y=(a ± σa)+(b ± σb)*X"], marker=:rect)
    plot!([x1-dX],[Y[1]], lc=:white, label=" ")   # adds space in legend
    plot!(X, st.a .+ st.b*X, color=:lightblue, ribbon=(σx95,σx95), label=false)
    plot!(xx, st.a .+ st.b*xx, color=:reds, ribbon=(σp,σm), label=false)

    #TODO: 95%-confidence interval computation for X-Y errors
    if all(≈(0), σX) && all(≈(0), σY)
        plot!(xx, st.a .+ st.b*xx, color=:blue, lw=0.5, xlabel="X", ylabel="Y", label="LinearFitXY")
        plot!(xlims=extrema(X))
    end

    scatter!(X, Y, msw=0.1, ms=1., msc=:lightgrey, xerror= σX, yerror= σY, label=false)
    scatter!(X, Y, msw=0.1, ms=1.5, mc=:blue, label=false)
    
    if !all(≈(0), σX) || !all(≈(0), σY)
        plot_covariance_ellipses!(X, Y, σX.^2, covXY, σY.^2; lc=:grey)
    end
end
