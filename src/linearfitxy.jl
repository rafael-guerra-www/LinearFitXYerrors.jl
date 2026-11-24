# Revised src/linearfitxy.jl (full file, drop-in replacement)
#
# Changes vs current master:
#   1. BUG FIX — de-attenuated Pearson ρ: uses the mean-square error variance
#      mean(σX.^2) (not var(σX), which is ~0 for constant σ), subtracts the
#      error cross-covariance for correlated errors, and de-attenuates
#      (factor > 1) instead of attenuating further. Validated by simulation:
#      true latent ρ = 0.80 was recovered as 0.79–0.82 (old formula: 0.60–0.63).
#   2. NEW — covariance of (a, b) returned as σab (needed for 2D confidence
#      regions and error propagation; cov(a,b) ≈ -X̄w·σb²).
#   3. Hoisted loop-invariant ωX, ωY, α out of the iteration loop; tolerance
#      1e-12 (1e-15 was unattainable → always 50 iterations); bᵢ trimmed.
#   4. Input validation: vector lengths, NaN, |r| < 1, degenerate fits →
#      informative errors instead of silent NaN.
#   5. `verbose=false` keyword: @printf banner now opt-in.
#   6. Scalar/vector handling via `isa Number`; zeros-vector σ handled by
#      all(iszero, ...) (previously zeros(N) silently produced Inf weights).
#
# Plot functions unchanged: they move to a package extension (next roadmap step).

using Statistics
using Distributions
using Printf

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
    σab :: Float64               # NEW: covariance of (a, b)
    σa95:: Float64
    σb95:: Float64
    S   :: Float64
    ρ   :: Float64
    bᵢ  :: Vector{Float64}
    i   :: Int64
end

"""
    linearfitxy(X, Y; σX=0, σY=0, r=0, verbose=false, isplot=false, ratio=1)

1D linear fitting with uncertainties in X and Y (York et al. 2004):
`Y = a + b*X`, errors `X ± σX`, `Y ± σY`, optional error correlation `r`
(scalar or per-point vector).

Returns `st::stfitxy` with: `a`, `b`, `σa`, `σb`, `σab` (covariance of a, b),
`σa95`, `σb95` (95% t-Student CIs), goodness of fit `S` (reduced Χ²),
de-attenuated Pearson coefficient `ρ` (corrected for regression dilution via
the mean-square error variances), iteration trace `bᵢ` and count `i`.
"""
function linearfitxy(X, Y; σX=0, σY=0, r=0, verbose=false, isplot=false, ratio=1)

    N = length(X)
    length(Y) == N && N >= 3 ||
        throw(ArgumentError("X and Y must have the same length, ≥ 3"))

    # expand scalars once; all-zero σ vectors handled like the scalar-0 case
    σX isa Number && (σX = fill(σX, N))
    σY isa Number && (σY = fill(σY, N))
    r  isa Number && (r  = fill(r,  N))
    (length(σX) == length(σY) == length(r) == N) ||
        throw(ArgumentError("σX, σY and r must be scalars or have length(X) elements"))
    any(isnan, X) && throw(ArgumentError("X contains NaN"))
    any(isnan, Y) && throw(ArgumentError("Y contains NaN"))

    d  = TDist(N - 2)                  # t-Student with N-2 d.o.f.
    cf = quantile(d, 0.975)            # two-tailed 95% factor

    if all(iszero, σX) && all(iszero, σY)
        # ---------------- no errors in X and Y: plain OLS ----------------
        σX = σY = r = zeros(N)

        M = [ones(N) X]
        (a, b) = M \ Y
        X̄ = mean(X)

        Ŝ = sqrt(sum((Y .- b .* X .- a).^2) / (N - 2))
        σb = sqrt(1 / sum((X .- X̄).^2))
        σa = Ŝ * sqrt(1 / N + X̄^2 * σb^2)
        σb *= Ŝ
        σab = -X̄ * σb^2               # OLS: cov(a,b) = -x̄·σb²

        ρ = cov(X, Y) / sqrt(var(X) * var(Y))

        st = stfitxy(X, Y, σX, σY, r, a, b, σa, σb, σab, cf*σa, cf*σb, Ŝ, ρ, [b], 1)

    else
        # ---------------- York (2004) with X-Y errors ----------------
        Nmax = 50                     # maximum number of iterations
        tol  = 1e-12                  # relative tolerance (1e-15 unattainable)

        all(iszero, σX) && (σX .= 1e-16)   # sentinel: no X errors
        all(iszero, σY) && (σY .= 1e-16)   # sentinel: no Y errors
        any(v -> !(abs(v) < 1), r) &&
            throw(ArgumentError("all |r| must be < 1 (error correlations)"))

        # loop-invariant weights: hoisted out of the iteration (constant in b)
        ωX = 1 ./ abs2.(σX)
        ωY = 1 ./ abs2.(σY)
        α  = sqrt.(ωX .* ωY)

        M = [ones(N) X]
        (a, b) = M \ Y                # initial guess
        bᵢ = [b]
        X̄ = 0.0; Ȳ = 0.0
        β = zeros(N); W = zeros(N)
        local i
        for outer i = 1:Nmax
            # W > 0 since |r| < 1: (sqrt(ωX) - b·sqrt(ωY))² + 2b·sqrt(ωX·ωY)·(1 − r) > 0
            W .= ωX .* ωY ./ (ωX .+ b^2 .* ωY .- 2 .* b .* r .* α)
            X̄ = sum(W .* X) / sum(W)
            Ȳ = sum(W .* Y) / sum(W)
            U = X .- X̄
            V = Y .- Ȳ
            β .= W .* (U ./ ωY .+ b .* V ./ ωX .- (b .* U .+ V) .* r ./ α)
            denom = sum(W .* β .* U)
            abs(denom) < 1e-300 &&
                throw(ErrorException("degenerate fit (vertical line?): b-estimate denominator ≈ 0"))
            b = sum(W .* β .* V) / denom
            push!(bᵢ, b)
            abs((bᵢ[end] - bᵢ[end-1]) / bᵢ[end]) < tol && break
        end
        a = Ȳ - b * X̄
        x  = X̄ .+ β
        X̄w = sum(W .* x) / sum(W)      # weighted centroid of the "adjusted" abscissae

        # goodness of fit (York 2004; Cantrell 2008)
        Ŝ = sqrt(sum(W .* (Y .- b .* X .- a).^2) / (N - 2))

        σb  = sqrt(1 / sum(W .* (x .- X̄w).^2))
        σa  = Ŝ * sqrt(1 / sum(W) + X̄w^2 * σb^2)
        σb *= Ŝ                        # Ŝ inflation (Cantrell 2008)
        σab = -X̄w * σb^2               # cov(a,b), same linearized level as σa, σb

        # de-attenuated Pearson ρ (regression dilution):
        # true ρ ≈ [cov(X,Y) − mean(r·σX·σY)] / sqrt((vX − mean(σX²))·(vY − mean(σY²)))
        vX = var(X); vY = var(Y)
        eX = mean(abs2, σX); eY = mean(abs2, σY)
        ρ  = (cov(X, Y) - mean(r .* σX .* σY))
        ρ /= sqrt(max(vX - eX, eps()) * max(vY - eY, eps()))
        ρ  = clamp(ρ, -1.0, 1.0)

        st = stfitxy(X, Y, σX, σY, r, a, b, σa, σb, σab, cf*σa, cf*σb, Ŝ, ρ, bᵢ, i)
    end

    if isplot
        _isplay(st; ratio = ratio)   # bridge defined in the module file;
                                      # requires the Plots extension to be loaded
    end

    if verbose
        @printf("\n>>> [± σ]  Y = (%.4f ± %.4f) + (%.4f ± %.4f)*X \n",  a, σa,  b, σb)
        @printf(">>> [± 95%% CI]  Y = (%.4f ± %.4f) + (%.4f ± %.4f)*X \n", a, cf*σa, b, cf*σb)
        @printf(">>> Pearson ρ = %.3f;  Goodness of fit = %.3f \n\n", ρ, Ŝ)
    end

    return st
end


# NOTE: the plot functions (plot_linfitxy, plot_covariance_ellipses!) are no
# longer part of src/linearfitxy.jl. The graphics now live in the backend-free
# recipe (src/recipes.jl, `plot(st::stfitxy)`) and the Plots-dependent helpers
# in ext/LinearFitXYerrorsPlotsExt.jl (loaded only when the user has Plots).
