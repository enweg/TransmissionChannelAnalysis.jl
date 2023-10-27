
@testset "Q constructors" begin
    q1 = TransmissionMechanisms.Q("x1")
    q2 = TransmissionMechanisms.Q("x2", -1.5)
    q3 = TransmissionMechanisms.Q(["x1", "x2"])
    q4 = TransmissionMechanisms.Q(["x1", "x2"], [-1, 2])

    @test q1.vars[1] == "x1"
    @test q2.vars[1] == "x2"
    @test q2.multiplier[1] == -1.5
    @test length(q3.vars) == 2
    @test all(q4.multiplier .== [-1, 2])
end

@testset "collect_terms" begin
    q = TransmissionMechanisms.Q(["x1", "x1"], [1, 1])  
    q = TransmissionMechanisms.collect_terms(q)
    @test length(q.vars) == 1
    @test q.vars[1] == "x1"
    @test q.multiplier[1] == 2

    q = TransmissionMechanisms.Q(["x1", "", "x1"], [1, 1, -1])  
    q = TransmissionMechanisms.collect_terms(q)
    @test length(q.vars) == 1
    @test q.vars[1] == ""
    @test q.multiplier[1] == 1
end

@testset "string_and" begin
    s = TransmissionMechanisms.string_and("", "x1")
    @test s == "x1"

    s = TransmissionMechanisms.string_and("x1", "x1")
    @test s == "x1"

    s = TransmissionMechanisms.string_and("x1 & x2", "x1")
    @test s == "x2 & x1"

    s = TransmissionMechanisms.string_and("!x1 & x2", "x1")
    @test s == "x2 & x1 & !x1"

    s = TransmissionMechanisms.string_and("x1", "x2")
    @test s == "x2 & x1"
end

@testset "AND" begin
    x = [TransmissionMechanisms.Q("x$i") for i = 1:5]

    q = x[1] & x[2]
    @test q.vars[1] == "x2 & x1"

    q = x[1] & TransmissionMechanisms.Q(["x1", "x2"], [1, 1])
    @test q.vars[1] == "x1"
    @test q.vars[2] == "x2 & x1"
    @test all(q.multiplier .== [1, 1])


    q = x[1] & TransmissionMechanisms.Q(["x1", "x2"], [1, -2])
    @test q.vars[1] == "x1"
    @test q.vars[2] == "x2 & x1"
    @test all(q.multiplier .== [1, -2])
end

@testset "OR" begin
    x = [TransmissionMechanisms.Q(i) for i = 1:5]

    q = x[1] | x[2]
    @test length(q.vars) == 3
    @test "x1" in q.vars
    @test "x2" in q.vars
    @test "x2 & x1" in q.vars

    @test q.multiplier[findfirst(==("x1"), q.vars)] == 1
    @test q.multiplier[findfirst(==("x2"), q.vars)] == 1
    @test q.multiplier[findfirst(==("x2 & x1"), q.vars)] == -1

    q = (x[1] & x[2]) | x[1]
    @test length(q.vars) == 1
    @test q.vars[1] == "x1"
    @test q.multiplier[1] == 1
end

@testset "NOT" begin
    x = [TransmissionMechanisms.Q(i) for i = 1:5]

    q = !x[1]
    @test q.vars[1] == "!x1"
    
    q = !(x[1] & x[2])
    @test length(q.vars) == 2
    @test "" in q.vars
    @test "x2 & x1" in q.vars
    @test q.multiplier[findfirst(==(""), q.vars)] == 1
    @test q.multiplier[findfirst(==("x2 & x1"), q.vars)] == -1

    q = !(x[1] | x[2])
    @test length(q.vars) == 4
    @test "" in q.vars
    @test "x1" in q.vars
    @test "x2" in q.vars
    @test "x2 & x1" in q.vars
    @test q.multiplier[findfirst(==(""), q.vars)] == 1
    @test q.multiplier[findfirst(==("x1"), q.vars)] == -1
    @test q.multiplier[findfirst(==("x2"), q.vars)] == -1
    @test q.multiplier[findfirst(==("x2 & x1"), q.vars)] == 1

    # An equivalent way of writing the above is Note that this representation
    # results in a shorter calculation. TODO: Can we make sure that the first
    # expression above also stops simplifying at !x1 & !x2?
    q = !x[1] & !x[2]
    @test length(q.vars) == 1
    @test "!x2 & !x1" in q.vars
    @test q.multiplier[1] == 1
end