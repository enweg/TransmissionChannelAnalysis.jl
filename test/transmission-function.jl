@testset "create_transmission_function" begin

    irfs = read_json_to_matrix("./simulated-svar-k3-p1/irfs.json")
    irfs_ortho = read_json_to_matrix("./simulated-svar-k3-p1/irfs_ortho.json")

    B = read_json_to_matrix("./simulated-svar-k3-p1/B.json")
    Qbb = read_json_to_matrix("./simulated-svar-k3-p1/Qbb.json")

    mat = test1(irfs, irfs_ortho, B, Qbb)
    @test isapprox(mat[:, 1], mat[:, 4]; atol = sqrt(eps()))
    mat = mat[3:end, :]
    @test isapprox(mat[:, 1], mat[:, 3]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 2]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 4]; atol = sqrt(eps()))

    mat = test2(irfs, irfs_ortho, B, Qbb)
    @test isapprox(mat[:, 1], mat[:, 4]; atol = sqrt(eps()))
    mat = mat[4:end, :]    
    @test isapprox(mat[:, 1], mat[:, 3]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 2]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 4]; atol = sqrt(eps()))

    mat = test3(irfs, irfs_ortho, B, Qbb)
    @test isapprox(mat[:, 1], mat[:, 4]; atol = sqrt(eps()))
    mat = mat[3:end, :]
    @test isapprox(mat[:, 1], mat[:, 3]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 2]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 4]; atol = sqrt(eps()))

    mat = test4(irfs, irfs_ortho, B, Qbb)
    @test isapprox(mat[:, 1], mat[:, 4]; atol = sqrt(eps()))
    mat = mat[4:end, :]
    @test isapprox(mat[:, 1], mat[:, 3]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 2]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 4]; atol = sqrt(eps()))

    mat = test5(irfs, irfs_ortho, B, Qbb)
    @test isapprox(mat[:, 1], mat[:, 4]; atol = sqrt(eps()))
    mat = mat[6:end, :]
    @test isapprox(mat[:, 1], mat[:, 3]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 2]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 4]; atol = sqrt(eps()))

    mat = test6(irfs, irfs_ortho, B, Qbb)
    @test isapprox(mat[:, 1], mat[:, 4]; atol = sqrt(eps()))
    mat = mat[6:end, :]
    @test isapprox(mat[:, 1], mat[:, 3]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 4]; atol = sqrt(eps()))

    mat = test7(irfs, irfs_ortho, B, Qbb)
    @test isapprox(mat[:, 1], mat[:, 3]; atol = sqrt(eps()))
    @test isapprox(mat[:, 1], mat[:, 4]; atol = sqrt(eps()))
end
