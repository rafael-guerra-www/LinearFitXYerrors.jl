# src/recipes.jl — backend-free figure geometry + Plots recipe
#
# Architecture (chosen so that the future Makie extension shares everything):
#   1. `plotdata(st)` computes ALL geometry once, as pure data (vectors,
#      NaN-separated segment lists, ellipse polygons). No backend, no Plots.
#   2. The Plots `@recipe` is a THIN renderer of that geometry.
#   3. The future Makie extension will be another thin renderer of the same
#      `plotdata` output (Band, errorbar-like segments via `lines!`, poly).
#
# Error bars are drawn as explicit NaN-separated `:path` segments instead of
# `xerror`/`yerror`: (a) error-bar attributes inside recipes are a known Plots
# problem area, (b) NaN-separated paths render identically in Plots and Makie.
#
# Fixed along the way: `maximum(...; dims=2)` returns an n×1 Matrix; without
# `vec()` the subsequent `.-` against an n-vector broadcast to an n×n Matrix,
# which made Plots draw garbage "vertical bars" (the bug seen on 2026-11-24).

using RecipesBase
using Statistics
using Printf
using Distributions: TDist, quantile

"""
    plotdata(st::stfitxy; level = 0.95) -> NamedTuple

Backend-free geometry of the standard `LinearFitXYerrors` figure.

Fields:
- `X`, `Y`: data sorted by X; `σX`, `σY`, `r`: per-point errors (sorted)
- `xx`, `fit`: fit-line path at the same abscissae
- `σp`, `σm`: upper/lower half-widths of the pointwise confidence band of the
  fitted line at confidence `level` (default 95%, the standard convention):
  `ŷ ± t · sqrt(σa² + 2x·σab + x²·σb²)`, exact under the linearized covariance
  (uses the σab field; the minimum-bandwidth point is the weighted centroid)
- `vbars`, `hbars`: `(x, y)` NaN-separated vertical/horizontal error-bar segments
- `ellipses`: `(x, y)` NaN-separated covariance-error-ellipse outlines
  (empty unless BOTH σX and σY are non-zero)
- `x1`, `x2`: x-axis extension limits
"""
function plotdata(st::stfitxy; level = 0.95)
    I  = sortperm(st.X)
    X  = st.X[I];  Y  = st.Y[I]
    σX = st.σX[I]; σY = st.σY[I]; r = st.r[I]

    dX = diff([extrema(X)...])[1] / 7          # extend x-axis by 1/7 each side
    x1, x2 = (-dX, dX) .+ extrema(X)
    xx = [x1; X; x2]
    fit = st.a .+ st.b .* xx

    # ---- pointwise confidence band of the line (standard: level = 0.95) ----
    # Var(ŷ(x)) = σa² + 2x·σab + x²·σb² (exact for the linearized covariance;
    # σab = -X̄w·σb² places the band minimum at the weighted centroid X̄w)
    N    = length(X)
    zcf  = quantile(TDist(N - 2), 1 - (1 - level) / 2)
    σŷ   = sqrt.(max.(st.σa^2 .+ 2 .* xx .* st.σab .+ (xx .^ 2) .* st.σb^2, 0.0))
    σp   = zcf .* σŷ
    σm   = copy(σp)

    # ---- error bars as explicit NaN-separated segments ----
    vx = Float64[]; vy = Float64[]
    hx = Float64[]; hy = Float64[]
    for i in eachindex(X)
        if !iszero(σY[i])
            append!(vx, (X[i], X[i], NaN))
            append!(vy, (Y[i] + σY[i], Y[i] - σY[i], NaN))
        end
        if !iszero(σX[i])
            append!(hx, (X[i] + σX[i], X[i] - σX[i], NaN))
            append!(hy, (Y[i], Y[i], NaN))
        end
    end

    # ---- covariance error ellipses (both X and Y errors; σX==0 would give
    # degenerate vertical segments duplicating the Y error bars) ----
    ex = Float64[]; ey = Float64[]
    if !all(iszero, σX) && !all(iszero, σY)
        ct = cos.(range(0, 2π; length = 73))
        stt = sin.(range(0, 2π; length = 73))
        for i in eachindex(X)
            aᵢ = σX[i]^2;  bᵢ = σX[i] * σY[i] * r[i];  cᵢ = σY[i]^2
            rᵢ = sqrt(((aᵢ - cᵢ)/2)^2 + bᵢ^2)
            λ₁ = sqrt(abs((aᵢ + cᵢ)/2 + rᵢ))
            λ₂ = sqrt(abs((aᵢ + cᵢ)/2 - rᵢ))
            θ  = (bᵢ == 0 && aᵢ ≥ cᵢ) ? 0.0 :
                 (bᵢ == 0 && aᵢ <  cᵢ) ? π/2 : atan(λ₁^2 - aᵢ, bᵢ)
            append!(ex, X[i] .+ λ₁ * cos(θ) .* ct .- λ₂ * sin(θ) .* stt)
            append!(ey, Y[i] .+ λ₁ * sin(θ) .* ct .+ λ₂ * cos(θ) .* stt)
            push!(ex, NaN); push!(ey, NaN)
        end
    end

    return (X = X, Y = Y, σX = σX, σY = σY,
            xx = xx, fit = fit, σp = σp, σm = σm,
            vbars = (x = vx, y = vy), hbars = (x = hx, y = hy),
            ellipses = (x = ex, y = ey), x1 = x1, x2 = x2)
end

# ============================ Plots recipe ============================

@recipe function f(st::stfitxy)
    d = plotdata(st)

    # ---------------- title, axes, legacy typography ----------------
    str1 = @sprintf("[± σ]  Y = (%.4f ± %.4f) + (%.4f ± %.4f)*X", st.a, st.σa, st.b, st.σb)
    str2 = @sprintf("\n[± 95%% CI]  Y = (%.4f ± %.4f) + (%.4f ± %.4f)*X", st.a, st.σa95, st.b, st.σb95)
    str3 = @sprintf("\nPearson ρ = %.2f; Goodness of fit = %.2f", st.ρ, st.S)
    title --> str1 * str2 * str3
    xguide --> "X"
    yguide --> "Y"
    legend --> :outerbottomright
    titlefontsize --> 8
    legendfontsize --> 6
    tickfontsize  --> 6
    guidefontsize --> 7
    fontfamily    --> "Computer Modern"
    xlims --> (d.x1, d.x2)

    # ---------------- covariance ellipses (behind everything) ----------------
    if !isempty(d.ellipses.x)
        @series begin
            linecolor    := :grey
            linewidth    := 0.3
            label        := ""
            primary      := false
            seriestype   := :path
            d.ellipses.x, d.ellipses.y
        end
    end

    # ---------------- error bars (explicit NaN-separated segments) ----------------
    if !isempty(d.vbars.x)
        @series begin
            linecolor    := :lightgrey
            linewidth    := 0.3
            label        := ""
            primary      := false
            seriestype   := :path
            d.vbars.x, d.vbars.y
        end
    end
    if !isempty(d.hbars.x)
        @series begin
            linecolor    := :lightgrey
            linewidth    := 0.3
            label        := ""
            primary      := false
            seriestype   := :path
            d.hbars.x, d.hbars.y
        end
    end

    # ---------------- fit line + pointwise 95% confidence band ----------------
    @series begin
        ribbon       := (d.σp, d.σm)
        fillcolor    := :lightblue
        linecolor    := :blue
        linewidth    := 0.5
        label        := "LinearFitXY (95% CI)"
        seriestype   := :path
        d.xx, d.fit
    end

    # ---------------- data points (markers only, no xerror/yerror) ----------------
    @series begin
        seriestype   := :scatter
        markersize   := 1.5
        markercolor  := :blue
        markerstrokecolor := :lightgrey
        label        := ""
        d.X, d.Y
    end
end
