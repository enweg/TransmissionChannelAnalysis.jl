@testset "create_transmission_function" begin

    irfs = read_json_to_matrix("./simulated-svar-k3-p1/irfs.json")
    irfs_ortho = read_json_to_matrix("./simulated-svar-k3-p1/irfs_ortho.json")

    B = read_json_to_matrix("./simulated-svar-k3-p1/B.json")
    Omega = read_json_to_matrix("./simulated-svar-k3-p1/Omega.json")

    mat = test1(irfs, irfs_ortho, B, Omega)
    # NOT conditions are never masked
    @test !any(isnan, mat)
    @test isapprox(mat[:, 1], mat[:, 4]; atol = sqrt(eps()))
    @test isapprox(mat[:, 2], mat[:, 3]; atol = sqrt(eps()))
    mat = mat[3:end, :]
    @test isapprox(mat[:, 1], mat[:, 3]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 2]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 4]; atol = sqrt(eps()))

    mat = test2(irfs, irfs_ortho, B, Omega)
    # x1 to x3 are not interpretable since all paths must go through x3
    @test all(isnan, mat[1:3, 3:4])
    mat = mat[4:end, :]    
    @test isapprox(mat[:, 1], mat[:, 3]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 2]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 4]; atol = sqrt(eps()))

    mat = test3(irfs, irfs_ortho, B, Omega)
    if TransmissionChannelAnalysis.REMOVE_CONTRADICTIONS[]
        @test isapprox(mat[:, 1], mat[:, 4]; atol = sqrt(eps()))
    else
        # without removing contradictions x2 remains an AND variable
        @test all(isnan, mat[1:2, 3:4])
    end
    mat = mat[3:end, :]
    @test isapprox(mat[:, 1], mat[:, 3]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 2]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 4]; atol = sqrt(eps()))

    mat = test4(irfs, irfs_ortho, B, Omega)
    @test all(isnan, mat[1:3, 3:4])
    mat = mat[4:end, :]
    @test isapprox(mat[:, 1], mat[:, 3]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 2]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 4]; atol = sqrt(eps()))

    mat = test5(irfs, irfs_ortho, B, Omega)
    @test all(isnan, mat[1:5, 3:4])
    mat = mat[6:end, :]
    @test isapprox(mat[:, 1], mat[:, 3]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 2]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 4]; atol = sqrt(eps()))

    mat = test6(irfs, irfs_ortho, B, Omega)
    @test all(isnan, mat[1:2, 3:4])
    @test isapprox(mat[3:end, 1], mat[3:end, 4]; atol = sqrt(eps()))
    mat = mat[6:end, :]
    @test isapprox(mat[:, 1], mat[:, 3]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 4]; atol = sqrt(eps()))

    mat = test7(irfs, irfs_ortho, B, Omega)
    @test isapprox(mat[:, 1], mat[:, 3]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 4]; atol = sqrt(eps()))
end
