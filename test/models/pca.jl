@testset "PCA" begin
    # 1. Factor loadings are orthonormal
    X = randn(10_000, 4)
    pca = PCA(X, 2)
    @test maximum(abs, pca.lambda' * pca.lambda - I) < 1e-6
    # 2. Factors are uncorrelated
    @test maximum(abs, tril(pca.F' * pca.F, -1)) < 1e-6

    # 3. If number of factors is the same as number of variables,
    # the same reconstructed X should be approximately the true X.
    X = randn(10_000, 4)
    pca = PCA(X, 4)
    Xhat = pca.F * pca.lambda'
    @test maximum(abs, Xhat - X) < sqrt(eps())

    # 4. Recovery of factors and loadings up to rotation
    # If Fhat are the recovered factors, then the rotation can be
    # recovered via H = (F' * F) \ (F' * Fhat). We can then test
    # whether we recover the correct loading.
    F = randn(10_000, 2)
    lambda = randn(4, 2)
    X = F * lambda'
    pca = PCA(X, 2)
    Fhat = pca.F
    H = (F' * F) \ (F' * Fhat)
    @test maximum(abs, H * pca.lambda' - lambda') < sqrt(eps())
end

@testset "PCA IC" begin
    Random.seed!(6150533)
    T = 10_000
    m = 3
    n = 20

    F = randn(T, m)
    lambda = randn(n, m)
    X = F * lambda'

    pca = PCA(X, m)
    @test TransmissionChannelAnalysis._msr(pca) < sqrt(eps())
    @test TransmissionChannelAnalysis._msr(pca, m) < sqrt(eps())
    TransmissionChannelAnalysis._msr(pca, 2)
    TransmissionChannelAnalysis._msr(pca, 1)

    BaiNgPC1(pca, 1, m)
    BaiNgPC2(pca, 1, m)
    BaiNgPC3(pca, 1, m)

    BaiNgIC1(pca, 1, m)
    BaiNgIC2(pca, 1, m)
    BaiNgIC3(pca, 1, m)

    pca = PCA(X, 1)
    pca, ic_table = select_factors!(pca, m + 2, BaiNgPC1)
    @test size(pca.F, 2) == m
    pca, ic_table = select_factors!(pca, m + 2, BaiNgPC2)
    @test size(pca.F, 2) == m
    pca, ic_table = select_factors!(pca, m + 2, BaiNgPC3)
    @test size(pca.F, 2) == m
    pca, ic_table = select_factors!(pca, m + 2, BaiNgIC1)
    @test size(pca.F, 2) == m
    pca, ic_table = select_factors!(pca, m + 2, BaiNgIC2)
    @test size(pca.F, 2) == m
    pca, ic_table = select_factors!(pca, m + 2, BaiNgIC3)
    @test size(pca.F, 2) == m

    # adding a tiny bit of errors
    Random.seed!(6150533)
    T = 10_000
    m = 3
    n = 20
    F = randn(T, m)
    lambda = randn(n, m)
    X = F * lambda' + 0.1 * randn(T, n)
    pca = PCA(X, 1)
    pca, ic_table = select_factors!(pca, m + 2, BaiNgPC1)
    @test size(pca.F, 2) == m
    pca, ic_table = select_factors!(pca, m + 2, BaiNgPC2)
    @test size(pca.F, 2) == m
    pca, ic_table = select_factors!(pca, m + 2, BaiNgPC3)
    @test size(pca.F, 2) == m
    pca, ic_table = select_factors!(pca, m + 2, BaiNgIC1)
    @test size(pca.F, 2) == m
    pca, ic_table = select_factors!(pca, m + 2, BaiNgIC2)
    @test size(pca.F, 2) == m
    pca, ic_table = select_factors!(pca, m + 2, BaiNgIC3)
    @test size(pca.F, 2) == m

    # Below is just a check if it also runs if there is some noise
    # No official test is being run.
    # Consider this an implementation test.
    X = F * lambda' + 0.1 * randn(T, n)
    pca = PCA(X, 1)
    pca, ic_table = select_factors!(pca, m + 2, BaiNgPC1)
end
