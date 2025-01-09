@testset "make_B make_Omega ordering=1:3" begin
    # structural model
    A0 = randn(3, 3)
    A0inv = inv(A0)
    As = [randn(3, 3) for _ in 1:2]
    Psis = [randn(3, 3) for _ in 1:2]
    Sigma = A0inv * A0inv'

    L, D = TransmissionChannelAnalysis.make_L_D(Sigma)
    O = zeros(3, 3)
    Q = A0 * inv(L)

    B = [
        (I - D*L) O O O;
        D*Q'*As[1] (I - D*L) O O;
        D*Q'*As[2] D*Q'*As[1] (I - D*L) O;
        O D*Q'*As[2] D*Q'*As[1] (I - D*L);
    ]

    Omega = [
        D*Q' O O O;
        D*Q'*Psis[1] D*Q' O O;
        D*Q'*Psis[2] D*Q'*Psis[1] D*Q' O;
        O D*Q'*Psis[2] D*Q'*Psis[1] D*Q'
    ]

    # reduced form model
    As = [A0inv * A for A in As]
    Psis = [A0inv * Psi * A0 for Psi in Psis]

    # testing
    B_test = TransmissionChannelAnalysis.make_B(As, Sigma, 1:3, 3)
    @test maximum(abs, B - B_test) < sqrt(eps())
    B_test = TransmissionChannelAnalysis.make_B(As, Sigma, 1:3, 1)
    @test maximum(abs, B_test - B[1:6, 1:6]) < sqrt(eps())

    Omega_test = TransmissionChannelAnalysis.make_Omega(A0inv, Psis, Sigma, 1:3, 3)
    @test maximum(abs, Omega_test - Omega) < sqrt(eps())
    Omega_test = TransmissionChannelAnalysis.make_Omega(A0inv, Psis, Sigma, 1:3, 1)
    @test maximum(abs, Omega_test - Omega[1:6, 1:6]) < sqrt(eps())
end


@testset "make_B make_Omega ordering=[3,1,2]" begin
    order = [3, 1, 2]
    T = TransmissionChannelAnalysis.permmatrix(order)
    # structural model
    A0 = randn(3, 3)
    A0inv = inv(A0)
    As = [randn(3, 3) for _ in 1:2]
    Psis = [randn(3, 3) for _ in 1:2]

    Sigma = A0inv * A0inv'

    L, D = TransmissionChannelAnalysis.make_L_D(T * Sigma *T')
    O = zeros(3, 3)
    Q = A0 * T' * inv(L)

    B = [
        (I - D*L) O O O;
        D*Q'*As[1]*T' (I - D*L) O O;
        D*Q'*As[2]*T' D*Q'*As[1]*T' (I - D*L) O;
        O D*Q'*As[2]*T' D*Q'*As[1]*T' (I - D*L);
    ]

    Omega = [
        D*Q' O O O;
        D*Q'*Psis[1] D*Q' O O;
        D*Q'*Psis[2] D*Q'*Psis[1] D*Q' O;
        O D*Q'*Psis[2] D*Q'*Psis[1] D*Q'
    ]

    # reduced form model
    As = [A0inv * A for A in As]
    Psis = [A0inv * Psi * A0 for Psi in Psis]

    # testing
    B_test = TransmissionChannelAnalysis.make_B(As, Sigma, order, 3)
    @test maximum(abs, B - B_test) < sqrt(eps())
    B_test = TransmissionChannelAnalysis.make_B(As, Sigma, order, 1)
    @test maximum(abs, B_test - B[1:6, 1:6]) < sqrt(eps())

    Omega_test = TransmissionChannelAnalysis.make_Omega(A0inv, Psis, Sigma, order, 3)
    @test maximum(abs, Omega_test - Omega) < sqrt(eps())
    Omega_test = TransmissionChannelAnalysis.make_Omega(A0inv, Psis, Sigma, order, 1)
    @test maximum(abs, Omega_test - Omega[1:6, 1:6]) < sqrt(eps())
end

@testset "make_B make_Omega ordering=[3,2,1]" begin
    order = 3:-1:1
    T = TransmissionChannelAnalysis.permmatrix(order)
    # structural model
    A0 = randn(3, 3)
    A0inv = inv(A0)
    As = [randn(3, 3) for _ in 1:2]
    Psis = [randn(3, 3) for _ in 1:2]

    Sigma = A0inv * A0inv'

    L, D = TransmissionChannelAnalysis.make_L_D(T * Sigma *T')
    O = zeros(3, 3)
    Q = A0 * T' * inv(L)

    B = [
        (I - D*L) O O O;
        D*Q'*As[1]*T' (I - D*L) O O;
        D*Q'*As[2]*T' D*Q'*As[1]*T' (I - D*L) O;
        O D*Q'*As[2]*T' D*Q'*As[1]*T' (I - D*L);
    ]

    Omega = [
        D*Q' O O O;
        D*Q'*Psis[1] D*Q' O O;
        D*Q'*Psis[2] D*Q'*Psis[1] D*Q' O;
        O D*Q'*Psis[2] D*Q'*Psis[1] D*Q'
    ]

    # reduced form model
    As = [A0inv * A for A in As]
    Psis = [A0inv * Psi * A0 for Psi in Psis]

    # testing
    B_test = TransmissionChannelAnalysis.make_B(As, Sigma, order, 3)
    @test maximum(abs, B - B_test) < sqrt(eps())
    B_test = TransmissionChannelAnalysis.make_B(As, Sigma, order, 1)
    @test maximum(abs, B_test - B[1:6, 1:6]) < sqrt(eps())

    Omega_test = TransmissionChannelAnalysis.make_Omega(A0inv, Psis, Sigma, order, 3)
    @test maximum(abs, Omega_test - Omega) < sqrt(eps())
    Omega_test = TransmissionChannelAnalysis.make_Omega(A0inv, Psis, Sigma, order, 1)
    @test maximum(abs, Omega_test - Omega[1:6, 1:6]) < sqrt(eps())
end



