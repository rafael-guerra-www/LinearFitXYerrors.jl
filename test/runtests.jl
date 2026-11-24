# test/runtests.jl — LinearFitXYerrors v0.2.0
# Structure:
#   1. Reference tests: the four README datasets with published values.
#      a, b, σa, σb, σa95, σb95, S are UNCHANGED vs v0.1.0 (regression guard:
#      the extension/cleanup work must not move the numbers). Only the ρ
#      expectations changed (v0.2.0 semantics: observed Pearson coefficient).
#   2. σab tests: new covariance field, exact external formula for OLS cases.
#   3. API/regression tests: the traps found during the v0.2.0 rework
#      (Int defaults → InexactError, zeros-vector σ, r range, NaN, lengths).
#   4. plotdata geometry tests (no plotting backend needed).
#   5. Plots extension smoke test (requires Plots in the test target).

using LinearFitXYerrors
using Test
using Statistics
using Plots                       # test target = ["Test", "Plots"]

@testset "LinearFitXYerrors.jl" begin

    @testset "0a - no errors (pba.ucdavis.edu data)" begin
        X = [16.4, 17.2, 17.6, 18.0, 18.2, 18.5]
        Y = [2.67, 2.75, 2.99, 3.14, 3.88, 4.23]
        st = linearfitxy(X, Y)
        @test st.a    ≈ -9.4628187535736
        @test st.b    ≈  0.7217838765009
        @test st.σa   ≈  3.6596191146816
        @test st.σb   ≈  0.2071823670235
        @test st.σa95 ≈ 10.1607315778461
        @test st.σb95 ≈  0.5752304688056
        @test st.S    ≈  0.3537301433210
        @test st.ρ    ≈  0.8672492974278     # plain Pearson, unchanged vs v0.1
        @test st.σab  ≈ -mean(X) * st.σb^2        # NEW field: OLS cov(a,b) = -x̄·σb²
        @test st.σab  ≈ -0.75761801107637
    end

    @testset "0b - no errors (Altman & Gardner 1988, Table I)" begin
        X = [91., 104, 107, 107, 106, 100, 92, 92, 105, 108]
        Y = [9.8, 7.4, 7.9, 8.3, 8.3, 9.0, 9.7, 8.8, 7.6, 6.9]
        st = linearfitxy(X, Y)
        @test st.a    ≈  20.1891143911439
        @test st.b    ≈ -0.1167896678967
        @test st.σa   ≈  2.6713186714425
        @test st.σb   ≈  0.0263407280356
        @test st.σa95 ≈  6.1600719027946
        @test st.σb95 ≈  0.0607418277744
        @test st.S    ≈  0.5484946692667
        @test st.ρ    ≈ -0.8430654646904     # plain Pearson, unchanged vs v0.1
        @test st.σab  ≈ -mean(X) * st.σb^2
        @test st.σab  ≈ -0.07021599608868
    end

    @testset "1 - uncorrelated errors (York 1966 / Pearson 1901)" begin
        X = [0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4]
        Y = [5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5]
        σX = 1 ./ sqrt.([1000., 1000, 500, 800, 200, 80,  60, 20, 1.8, 1])
        σY = 1 ./ sqrt.([1., 1.8, 4, 8, 20, 20, 70, 70, 100, 500])
        r = 0 * X
        st = linearfitxy(X, Y; σX = σX, σY = σY, r = r)
        @test st.a    ≈  5.47991022403286
        @test st.b    ≈ -0.48053340744620
        @test st.σa   ≈  0.35924652255111
        @test st.σb   ≈  0.07062026952877
        @test st.σa95 ≈  0.82842396656058
        @test st.σb95 ≈  0.16285063356258
        @test st.S    ≈  1.21790564053940
        @test st.ρ    ≈ -0.97647522267457    # v0.2.0: observed Pearson coefficient
                                             # (v0.1.0 reported -0.94325300584483)
        @test st.σab  < 0                     # b>0 => weighted centroid > 0 scaling
        @test all(isfinite, [st.a, st.b, st.σa, st.σb, st.σab])
    end

    @testset "2 - correlated errors (Amen 2012, Table 3.1)" begin
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
        μₑ = M[:, 2];  η = M[:, 3]
        σμₑ = sqrt.(M[:, 4]);  ση = sqrt.(M[:, 5])
        cov_μₑ_η = M[:, 6]
        r = cov_μₑ_η ./ (σμₑ .* ση)
        st = linearfitxy(μₑ, η; σX = σμₑ, σY = ση, r = r, isplot = false)
        @test st.a  ≈ -21.9553309045508
        @test st.b  ≈   1.1762288806682
        @test st.σa ≈   1.2846101711080
        @test st.σb ≈   0.0578260394408
        @test st.S  ≈   3.4885539308643
        @test st.ρ  ≈   0.89217586033265    # v0.2.0: observed Pearson coefficient
                                            # (v0.1.0 reported 0.8435087435271)
    end

    @testset "API and v0.2.0 regression traps" begin
        # pure-Integer defaults (regression: InexactError on 1e-16 sentinel)
        st = linearfitxy([1, 2, 3], [2, 4, 6])
        @test st.b ≈ 2.0 atol = 1e-9        # Y = 2X exactly: b = 2, a = 0
        @test st.a ≈ 0.0 atol = 1e-9
        @test st.σab isa Float64

        # single-side errors: σX = 0 (scalar), σY vector
        X = collect(1.0:5.0)
        Y = 2.0 .+ 0.5 .* X .+ 0.1 .* randn(5)
        st = linearfitxy(X, Y; σY = 0.1 .* ones(5))
        @test all(isfinite, [st.a, st.b, st.σa, st.σb, st.σab, st.S, st.ρ])

        # σX given as a zeros VECTOR (previously produced Inf weights silently)
        st = linearfitxy(X, Y; σX = zeros(5), σY = 0.1 .* ones(5))
        @test all(isfinite, [st.a, st.b, st.σb, st.S])

        # input validation
        @test_throws ArgumentError linearfitxy([1., 2, 3], [1., 2])           # length mismatch
        @test_throws ArgumentError linearfitxy([1., 2, 3], [1., 2, 3];      # σX wrong length
                                                σX = [1., 2])
        @test_throws ArgumentError linearfitxy([1., 2, 3], [1., 2, NaN])     # NaN input
        @test_throws ArgumentError linearfitxy([1., 2, 3], [2., 4, 6];       # |r| ≥ 1
                                                σX = 0.1, σY = 0.1, r = 1.0)

        # verbose kwarg accepted, and stfitxy is exported
        st = linearfitxy([1., 2, 3], [2., 4, 6]; verbose = false)
        @test st isa LinearFitXYerrors.stfitxy
    end

    @testset "plotdata geometry (no plotting backend needed)" begin
        # Y-errors only
        X = collect(1.0:10.0)
        Y = 2.0 .+ 0.5 .* X .+ 0.1 .* randn(10)
        σY = 0.1 .* ones(10)
        st = linearfitxy(X, Y; σY = σY)
        d = LinearFitXYerrors.plotdata(st)

        # vertical error-bar segments: 3 entries per point (top, bottom, NaN)
        @test length(d.vbars.x) == 3 * 10
        @test count(isnan, d.vbars.y) == 10
        # each segment spans exactly ±σY around its point (alternating top/bottom)
        seg = filter(!isnan, d.vbars.y)
        @test seg[1:2:end] ≈ Y .+ σY atol = 1e-12
        @test seg[2:2:end] ≈ Y .- σY atol = 1e-12

        # no horizontal bars, no ellipses (σX = 0)
        @test isempty(d.hbars.x)
        @test isempty(d.ellipses.x)

        # confidence band: symmetric, positive, finite
        @test d.σp == d.σm
        @test all(>(0), d.σp)
        @test all(isfinite, d.σp)

        # both-side errors: ellipses emitted (Amen-style correlated errors)
        Xa = [22.0, 23.5, 21.0]; Ya = [5.0, 8.0, 3.0]
        st = linearfitxy(Xa, Ya; σX = [0.3, 0.5, 0.2], σY = [0.4, 0.6, 0.3], r = 0.5)
        d = LinearFitXYerrors.plotdata(st)
        @test !isempty(d.ellipses.x)
        @test !isempty(d.hbars.x)
    end

    @testset "Plots extension smoke test" begin
        X = collect(1.0:10.0)
        Y = 2.0 .+ 0.5 .* X .+ 0.1 .* randn(10)
        st = linearfitxy(X, Y; σY = 0.1 .* ones(10))
        p = plot(st)
        @test p isa Plots.Plot
        # σY-only case: exactly 3 series — CI band + fit line, vertical error
        # bars, scatter. (Ellipses & horizontal bars are correctly absent
        # when σX = 0; with both-side errors the count rises to 5.)
        @test p.n == 3
    end

end
