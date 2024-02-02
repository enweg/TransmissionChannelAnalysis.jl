function test1(irfs, irfs_ortho, B, Qbb)
    cond = make_condition("!x2")
    effect = transmission(1, B, Qbb, cond)
    effect_irfs = transmission(1, irfs, irfs_ortho, cond; method = :irfs)

    B_tilde = copy(B)
    Qbb_tilde = copy(Qbb)
    Qbb_tilde[2, 1] = 0
    B_tilde[2, :] .= 0

    manual_2 = (inv(I - B_tilde)*Qbb_tilde)[:, 1]
    manual_1 = irfs[:, 1] - irfs[2, 1] * irfs_ortho[:, 2] / irfs_ortho[2, 2];
    manual_2 = manual_2[:, 1]
    return hcat(manual_1, manual_2, effect, effect_irfs)
end

function test2(irfs, irfs_ortho, B, Qbb)
    manual_1 = irfs[3,1] .* irfs_ortho[:,3] ./ irfs_ortho[3, 3] - 
        irfs[2,1] .* irfs_ortho[3,2] ./ irfs_ortho[2,2] .* irfs_ortho[:,3] ./ irfs_ortho[3,3]
    
    B_tilde = copy(B)
    Qbb_tilde = copy(Qbb)
    Qbb_tilde[4, 1] = 0
    B_tilde[4, [1, 2, 4]] .= 0

    manual_2_part1 = inv(I - B_tilde) * Qbb_tilde
    
    B_tilde = copy(B)
    Qbb_tilde = copy(Qbb)
    Qbb_tilde[[3,4], 1] .= 0
    B_tilde[3, 1] = 0
    B_tilde[4, [1, 2, 4]] .= 0

    manual_2_part2 = inv(I - B_tilde) * Qbb_tilde
    manual_2 = manual_2_part1 - manual_2_part2
    manual_2 = manual_2[:, 1]
    
    cond = make_condition("x3 & !x2")
    effect = transmission(1, B, Qbb, cond)
    effect_irfs = transmission(1, irfs, irfs_ortho, cond; method = :irfs)

    
    return hcat(manual_1, manual_2, effect, effect_irfs)
end

function test3(irfs, irfs_ortho, B, Qbb)
    cond = make_condition("x2 & !x2")
    effect = transmission(1, B, Qbb, cond)
    effect_irfs = transmission(1, irfs, irfs_ortho, cond; method = :irfs)
    return hcat(zeros(size(irfs, 1), 2), effect, effect_irfs)
end

function test4(irfs, irfs_ortho, B, Qbb)
    manual_1 = irfs[2, 1] .* irfs_ortho[:, 2] ./ irfs_ortho[2, 2] + 
        irfs[3, 1] .* irfs_ortho[:, 3] ./ irfs_ortho[3, 3] - 
        2 * irfs[2, 1] .* irfs_ortho[3, 2] ./ irfs_ortho[2, 2] .* irfs_ortho[:, 3] ./ irfs_ortho[3, 3]
    
    B_tilde = copy(B)
    Qbb_tilde = copy(Qbb)
    Qbb_tilde[[3, 4], 1] .= 0
    B_tilde[3:end, 1] .= 0

    manual_2_part1 = inv(I - B_tilde) * Qbb_tilde 

    B_tilde = copy(B)
    Qbb_tilde = copy(Qbb)
    Qbb_tilde[4,1] = 0
    B_tilde[4:end, [1, 2, 4]] .= 0

    manual_2_part2 = inv(I - B_tilde) * Qbb_tilde

    B_tilde = copy(B)
    Qbb_tilde = copy(Qbb)
    Qbb_tilde[[3, 4], 1] .= 0
    B_tilde[[3, 4], 1] .= 0
    B_tilde[4:end, [1, 2]] .= 0

    manual_2_part3 = inv(I - B_tilde) * Qbb_tilde

    manual_2 = manual_2_part1 + manual_2_part2 - 2 * manual_2_part3
    manual_2 = manual_2[:, 1]
    
    cond = make_condition("((x2 & !x3) | (!x2 & x3))")
    effect = transmission(1, B, Qbb, cond)
    effect_irfs = transmission(1, irfs, irfs_ortho, cond; method = :irfs)
    
    return hcat(manual_1, manual_2, effect, effect_irfs)
end

function test5(irfs, irfs_ortho, B, Qbb)
    
    manual_1 = irfs[2, 1] .* irfs_ortho[:, 2] ./ irfs_ortho[2, 2] + 
        irfs[5, 1] .* irfs_ortho[:, 5] ./ irfs_ortho[5, 5] - 
        irfs[2, 1] .* irfs_ortho[5, 2] ./ irfs_ortho[2, 2] .* irfs_ortho[:, 5] ./ irfs_ortho[5, 5]
    
    B_tilde = copy(B)
    Qbb_tilde = copy(Qbb)
    Qbb_tilde[3:end, 1] .= 0
    B_tilde[3:end, 1] .= 0
    manual_2_part1 = inv(I - B_tilde) * Qbb_tilde
    
    B_tilde = copy(B)
    Qbb_tilde = copy(Qbb)
    Qbb_tilde[6:end, 1] .= 0
    B_tilde[6:end, 1:4] .= 0
    manual_2_part2 = inv(I - B_tilde) * Qbb_tilde
    
    B_tilde = copy(B)
    Qbb_tilde = copy(Qbb)
    Qbb_tilde[3:end, 1] .= 0
    B_tilde[3:end, 1] .= 0
    B_tilde[6:end, 1:4] .= 0
    manual_2_part3 = inv(I - B_tilde) * Qbb_tilde
    
    manual_2 = (manual_2_part1 + manual_2_part2 - manual_2_part3)[:, 1]
    
    cond = make_condition("x2 | x5")
    effect = transmission(1, B, Qbb, cond)
    effect_irfs = transmission(1, irfs, irfs_ortho, cond; method = :irfs)
    
    
    return hcat(manual_1, manual_2, effect, effect_irfs)
end

function test6(irfs, irfs_ortho, B, Qbb)
    
    manual_1 = irfs[2, 1] .* irfs_ortho[:, 2] ./ irfs_ortho[2, 2] - 
        irfs[2, 1] .* irfs_ortho[3, 2] ./ irfs_ortho[2, 2] .* irfs_ortho[:, 3] ./ irfs_ortho[3, 3] - 
        irfs[2, 1] .* irfs_ortho[4, 2] ./ irfs_ortho[2, 2] .* irfs_ortho[:, 4] ./ irfs_ortho[4, 4] + 
        irfs[2, 1] .* irfs_ortho[3, 2] ./ irfs_ortho[2, 2] .* irfs_ortho[4, 3] ./ irfs_ortho[3, 3] .* irfs_ortho[:, 4] ./ irfs_ortho[4, 4] - 
        irfs[2, 1] .* irfs_ortho[5, 2] ./ irfs_ortho[2, 2] .* irfs_ortho[:, 5] ./ irfs_ortho[5, 5] + 
        irfs[2, 1] .* irfs_ortho[3, 2] ./ irfs_ortho[2, 2] .* irfs_ortho[5, 3] ./ irfs_ortho[3, 3] .* irfs_ortho[:, 5] ./ irfs_ortho[5, 5] + 
        irfs[2, 1] .* irfs_ortho[4, 2] ./ irfs_ortho[2, 2] .* irfs_ortho[5, 4] ./ irfs_ortho[4, 4] .* irfs_ortho[:, 5] ./ irfs_ortho[5, 5] - 
        irfs[2, 1] .* irfs_ortho[3, 2] ./ irfs_ortho[2, 2] .* irfs_ortho[4, 3] ./ irfs_ortho[3, 3] .* irfs_ortho[5, 4] ./ irfs_ortho[4, 4] .* irfs_ortho[:, 5] ./ irfs_ortho[5, 5]
    
    # too many terms to calculate it using the second method
    manual_2 = fill(NaN, size(irfs, 1))
    
    cond = make_condition("x2 & !x3 & !x4 & !x5")
    effect = transmission(1, B, Qbb, cond)
    effect_irfs = transmission(1, irfs, irfs_ortho, cond; method = :irfs)
    
    return hcat(manual_1, manual_2, effect, effect_irfs)
end

function test7(irfs, irfs_ortho, B, Qbb)
    manual_1 = irfs[:, 1]
    manual_2 = irfs[:, 1]

    cond = make_condition("(x1 | x2 | x3) | !(x1 | x2 | x3)")
    effect = transmission(1, B, Qbb, cond)
    effect_irfs = transmission(1, irfs, irfs_ortho, cond; method = :irfs)

    return hcat(manual_1, manual_2, effect, effect_irfs)
end

function read_json_to_matrix(file_name)
    f = read(file_name, String)
    mat = reduce(hcat, JSON.parse(f))
    return Float64.(mat)
end

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
