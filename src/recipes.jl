# src/recipes.jl — plot recipe for `stfitxy`
# Defines what `plot(st)` draws with any Recipes-based backend (Plots.jl):
# fit line with (± σa, ± σb) ribbon, optional 95%-CI band (no-errors case),
# covariance-error ellipses (X-Y errors case), and the data points with
# their X-Y error bars. Visual output matches the legacy plot_linfitxy.
# RecipesBase is backend-free: nothing here requires Plots to be loaded.

using RecipesBase
using Statistics
using Printf
using Distributions: TDist, quantile

@recipe function f(st::stfitxy)
    # ---------------- title and axes ----------------
    str1 = @sprintf("[± σ]  Y = (%.4f ± %.4f) + (%.4f ± %.4f)*X", st.a, st.σa, st.b, st.σb)
    str2 = @sprintf("\n[± 95%% CI]  Y = (%.4f ± %.4f) + (%.4f ± %.4f)*X", st.a, st.σa95, st.b, st.σb95)
    str3 = @sprintf("\nPearson ρ = %.2f; Goodness of fit = %.2f", st.ρ, st.S)
    title --> str1 * str2 * str3
    xguide --> "X"
    yguide --> "Y"
    legend --> :outerbottomright

    # ---------------- data sorted by X ----------------
    I    = sortperm(st.X)
    X    = st.X[I];  Y    = st.Y[I]
    σX   = st.σX[I]; σY   = st.Y[I]
    noerr = all(≈(0), σX) && all(≈(0), σY)

    dX = diff([extrema(X)...])[1] / 7       # extend x-axis by 1/7 on each side
    x1, x2 = (-dX, dX) .+ extrema(X)
    xx = [x1; X; x2]
    xlims --> (x1, x2)

    # ---------------- fit line + (±σa, ±σb) ribbon ----------------
    tl, bl = (st.a - st.σa) .+ (st.b + st.σb) .* xx, (st.a + st.σa) .+ (st.b - st.σb) .* xx
    σp = maximum([tl bl]; dims=2) .- (st.a .+ st.b .* xx)
    σm = (st.a .+ st.b .* xx) .- minimum([tl bl]; dims=2)

    @series begin
        ribbon       := (σp, σm)
        fillcolor    := :lightblue
        linecolor    := :blue
        linewidth    := 0.5
        label        := "LinearFitXY"
        seriestype   := :path
        xx, st.a .+ st.b .* xx
    end

    # ---------------- 95% confidence band (no-errors case only) ----------------
    if noerr
        N  = length(X)
        cf = quantile(TDist(N - 2), 0.975)
        σx95 = cf .* st.S .* sqrt.(1/N .+ (X .- mean(X)).^2 ./ var(X) ./ (N - 1))
        @series begin
            ribbon       := (σx95, σx95)
            fillcolor    := :reds
            linecolor    := :reds
            label        := "95% confidence interval"
            seriestype   := :path
            X, st.a .+ st.b .* X
        end
    else
        # ---------------- covariance error ellipses (X-Y errors case) ----------------
        # eigen-decomposition of the 2×2 covariance (cookierobotics.com/007)
        covXY = (σX .* σY .* st.r[I])
        t  = range(0, 2π; length = 73)
        ct = cos.(t); stt = sin.(t)
        for i in eachindex(X)
            aᵢ = σX[i]^2;  bᵢ = covXY[i];  cᵢ = σY[i]^2
            rᵢ = sqrt(((aᵢ - cᵢ)/2)^2 + bᵢ^2)
            λ₁ = sqrt(abs((aᵢ + cᵢ)/2 + rᵢ))
            λ₂ = sqrt(abs((aᵢ + cᵢ)/2 - rᵢ))
            θ  = (bᵢ == 0 && aᵢ ≥ cᵢ) ? 0.0 :
                 (bᵢ == 0 && aᵢ <  cᵢ) ? π/2 : atan(λ₁^2 - aᵢ, bᵢ)
            xe = X[i] .+ λ₁ * cos(θ) .* ct .- λ₂ * sin(θ) .* stt
            ye = Y[i] .+ λ₁ * sin(θ) .* ct .+ λ₂ * cos(θ) .* stt
            @series begin
                linecolor    := :grey
                linewidth    := 0.3
                label        := ""
                primary      := false      # excluded from legend/limits logic
                seriestype   := :path
                xe, ye
            end
        end
    end

    # ---------------- data points with X-Y error bars ----------------
    @series begin
        seriestype   := :scatter
        xerror       := σX
        yerror       := σY
        markersize   := 1.5
        markercolor  := :blue
        markerstrokecolor := :lightgrey
        label        := ""
        X, Y
    end
end
