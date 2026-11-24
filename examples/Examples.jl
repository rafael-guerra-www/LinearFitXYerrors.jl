# Updated examples/ scripts for v0.2.0
# Conventions for all figures: use plot(st) (the recipe) rather than the legacy
# isplot=true; explicit size/dpi for crisp README rendering; savefig for upload.
# New vs v0.1.0 figures: 95% CI band everywhere, exact error bars, ellipses only
# for two-side errors, and examples 1/4 overlay the single-side fits to display
# the regression-dilution effect (b(Y-errors only) vs b(XY errors)).

##############   EXAMPLE-0  (examples/example0.jl)   ##############
# No errors in X nor in Y: the 95% CI band of the fitted line is drawn

using LinearFitXYerrors
using Plots; gr()

# Example-0a: http://pba.ucdavis.edu/files/45007.pdf
X = [16.4, 17.2, 17.6, 18.0, 18.2, 18.5]
Y = [2.67, 2.75, 2.99, 3.14, 3.88, 4.23]
stxy = linearfitxy(X, Y)
plot(stxy, size = (600, 400), dpi = 300)
savefig(raw"examples\Example0a_LinearFitXYerrors.png")

# Example-0b: Altman & Gardner [1988], Table I
X = [91., 104, 107, 107, 106, 100, 92, 92, 105, 108]
Y = [9.8, 7.4, 7.9, 8.3, 8.3, 9.0, 9.7, 8.8, 7.6, 6.9]
stxy = linearfitxy(X, Y)
plot(stxy, size = (600, 400), dpi = 300)
savefig(raw"examples\Example0b_LinearFitXYerrors.png")


##############   EXAMPLE-1  (examples/example1.jl)   ##############
# Uncorrelated errors in X and Y, PLUS the single-side-error overlays
# (Y errors only / X errors only): shows the regression-dilution effect.

using LinearFitXYerrors
using Plots; gr()

# Tables I and II from York [1966] (data from Pearson 1901; Cantrell 2008)
X = [0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4]
Y = [5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5]
σX = 1 ./ sqrt.([1000., 1000, 500, 800, 200, 80,  60, 20, 1.8, 1])
σY = 1 ./ sqrt.([1., 1.8, 4, 8, 20, 20, 70, 70, 100, 500])
r = 0 * X

# full XY-errors fit (ellipses + both-side bars)
stxy = linearfitxy(X, Y; σX = σX, σY = σY, r = r)

# single-side fits: turn OFF errors in X (Y errors only) or in Y (X errors only)
sty = linearfitxy(X, Y; σX = 0, σY = σY, r = r)     # Y errors only
stx = linearfitxy(X, Y; σX = σX, σY = 0, r = r)     # X errors only

p = plot(stxy, size = (650, 420), dpi = 300)
dX = diff([extrema(X)...])[1] / 7
x1, x2 = (-dX, dX) .+ extrema(X)
xx = [x1; X; x2]
plot!(p, xx, sty.a .+ sty.b .* xx, color = :lime,   lw = 0.5, label = "LinearFitXY (Y errors only)")
plot!(p, xx, stx.a .+ stx.b .* xx, color = :orange, lw = 0.5, label = "LinearFitXY (X errors only)")
savefig(p, raw"examples\Example1_LinearFitXYerrors.png")


##############   EXAMPLE-2  (examples/example2.jl)   ##############
# Correlated errors in X and Y (Amen [2012], Table 3.1; r ≈ -0.93 average)

using LinearFitXYerrors
using Plots; gr()

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

μₑ = M[:, 2];    η = M[:, 3]
σμₑ = sqrt.(M[:, 4]);    ση = sqrt.(M[:, 5])
r = M[:, 6] ./ (σμₑ .* ση)      # per-point error correlation (includes negatives)

stxy = linearfitxy(μₑ, η; σX = σμₑ, σY = ση, r = r)
plot(stxy, size = (650, 420), dpi = 300)
savefig(raw"examples\Example2_LinearFitXYerrors.png")


##############   EXAMPLE-3  (examples/example3.jl)   ##############
# Uncorrelated vs correlated errors (Mahon [1996], Table 3)

using LinearFitXYerrors
using Plots; gr()

M = [0.037   0.0080
     0.035   0.0084
     0.032   0.0100
     0.040   0.0085
     0.013   0.0270
     0.038   0.0071
     0.042   0.0043
     0.030   0.0160]

X = M[:, 1];    Y = M[:, 2]
σX = 0.03 .* X
σY = 0.10 .* Y
r1 = 0          # uncorrelated errors
r2 = 0.7071     # correlated errors

stxy1 = linearfitxy(X, Y; σX = σX, σY = σY, r = r1)
plot(stxy1, size = (650, 420), dpi = 300)
savefig(raw"examples\Example3_LinearFitXYerrors.png")

# correlated errors, same data: only the ellipses and the fit change
stxy2 = linearfitxy(X, Y; σX = σX, σY = σY, r = r2)
plot(stxy2, size = (650, 420), dpi = 300)
savefig(raw"examples\Example3b_LinearFitXYerrors.png")


##############   EXAMPLE-4  (examples/example4.jl)   ##############
# NEW in v0.2.0: single-side errors made explicit, compared side by side.
# Fixed synthetic dataset with KNOWN ground truth (a=2, b=1.5), so each figure
# shows how far the recovered line is from the true line when errors on one
# side are ignored vs declared:
#   - Fig 4a: Y errors only declared (no ellipses, vertical bars, CI band)
#   - Fig 4b: X errors only declared (horizontal bars)
#   - Fig 4c: both declared (ellipses) ; overlays of 4a/4b reveal the bias of
#     ignoring errors in X (attenuated |b|) or in Y.

using LinearFitXYerrors
using Plots; gr()
using Random

Random.seed!(42)
N  = 30
Xₜ = range(1, 10; length = N)               # true (error-free) X
Yₜ = 2.0 .+ 1.5 .* Xₜ                       # true line: a=2, b=1.5
X  = Xₜ .+ 0.4 .* randn(N)                  # observed X (with error)
Y  = Yₜ .+ 0.4 .* randn(N)                  # observed Y (with error)
σX = 0.4 .* ones(N)                          # correctly declared errors
σY = 0.4 .* ones(N)

sty  = linearfitxy(X, Y; σX = 0,   σY = σY)  # Y errors only
stx  = linearfitxy(X, Y; σX = σX, σY = 0)    # X errors only
stxy = linearfitxy(X, Y; σX = σX, σY = σY)   # both sides declared

dX = diff([extrema(X)...])[1] / 7
x1, x2 = (-dX, dX) .+ extrema(X)
xx = [x1; X; x2]
overlay(p, st; label) = plot!(p, xx, st.a .+ st.b .* xx, lw = 0.5, label = label)

p = plot(sty, size = (650, 420), dpi = 300,
         title = "Errors in Y only — attenuated slope b = $(round(sty.b, digits=3)) (true: 1.5)")
overlay(p, stxy, label = "Fit with XY errors declared")
savefig(p, raw"examples\Example4a_LinearFitXYerrors.png")

p = plot(stx, size = (650, 420), dpi = 300,
         title = "Errors in X only — attenuated slope b = $(round(stx.b, digits=3)) (true: 1.5)")
overlay(p, stxy, label = "Fit with XY errors declared")
savefig(p, raw"examples\Example4b_LinearFitXYerrors.png")

p = plot(stxy, size = (650, 420), dpi = 300,
         title = "XY errors declared — recovered b = $(round(stxy.b, digits=3)) (true: 1.5)")
overlay(p, sty, label = "Fit ignoring errors in X")
overlay(p, stx, label = "Fit ignoring errors in Y")
savefig(p, raw"examples\Example4c_LinearFitXYerrors.png")
