using LinearFitXYerrors
using Test

@testset "LinearFitXYerrors.jl" begin


    # 0a - TEST DATA WITH NO ERRORS
    # INPUT DATA:  http://pba.ucdavis.edu/files/45007.pdf

    X = [16.4, 17.2, 17.6, 18.0, 18.2, 18.5]
    Y = [2.67, 2.75, 2.99, 3.14, 3.88, 4.23]

    stxy = linearfitxy(X,Y)

    @test   stxy.a    ≈ -9.4628187535736
    @test   stxy.b    ≈  0.7217838765009
    @test   stxy.σa   ≈  3.6596191146816
    @test   stxy.σb   ≈  0.2071823670235
    @test   stxy.σa95 ≈ 10.1607315778461
    @test   stxy.σb95 ≈  0.5752304688056
    @test   stxy.S    ≈  0.3537301433210
    @test   stxy.ρ    ≈  0.8672492974278


    # 0b - TEST DATA WITH NO ERRORS
    # INPUT DATA:  Table I from Altman and Gardner [1988]

    X = [91., 104, 107, 107, 106, 100, 92, 92, 105, 108]
    Y = [9.8, 7.4, 7.9, 8.3, 8.3, 9.0, 9.7, 8.8, 7.6, 6.9]

    stxy = linearfitxy(X,Y)
  
    @test   stxy.a    ≈ 20.1891143911439
    @test   stxy.b    ≈ -0.1167896678967
    @test   stxy.σa   ≈  2.6713186714425
    @test   stxy.σb   ≈  0.0263407280356
    @test   stxy.σa95 ≈  6.1600719027946
    @test   stxy.σb95 ≈  0.0607418277744
    @test   stxy.S    ≈  0.5484946692667
    @test   stxy.ρ    ≈ -0.8430654646904



    # 1 - TEST UNCORRELATED ERRORS
    # INPUT DATA:
    # Tables I and II from York [1966], data taken from Pearson (1901)
    # also in Cantrell (2008)

    X = [0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4]
    Y = [5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5]

    # the weights ωX, ωY are the inverse of the variances
    σX = 1 ./ sqrt.([1000., 1000, 500, 800, 200, 80,  60, 20, 1.8, 1])
    σY = 1 ./ sqrt.([1., 1.8, 4, 8, 20, 20, 70, 70, 100, 500])
    r = 0*X;

    stxy = linearfitxy(X,Y; σX=σX, σY=σY, r=r)

    @test   stxy.a  ≈  5.47991022403286
    @test   stxy.b  ≈ -0.48053340744620
    @test   stxy.σa ≈ 0.359246522551112
    @test   stxy.σb ≈ 0.070620269528771
    @test   stxy.S  ≈ 1.217905640539398
    @test   stxy.ρ  ≈ -0.94325300584483



    # 2 - TEST CORRELATED ERRORS

    # INPUT DATA: Table 3.1 from Amen, S. [2012]

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


    # COMPUTE:
    stxy = linearfitxy(μₑ, η; σX=σμₑ, σY=ση, r=r, isplot=false)

    @test   stxy.a  ≈ -21.9553309045508
    @test   stxy.b  ≈   1.1762288806682
    @test   stxy.σa ≈   1.2846101711080
    @test   stxy.σb ≈   0.0578260394408
    @test   stxy.S  ≈   3.4885539308643
    @test   stxy.ρ  ≈   0.8435087435271

end
