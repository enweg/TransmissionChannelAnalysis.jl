@testset "collect_and_terms" begin
    s = "x2 & x3 & !x2 & x5"
    cond = make_condition(s)
    and_terms = TransmissionMechanisms.collect_and_terms(cond)
    @test all(string.(and_terms) .== ["x5", "x3", "x2", "!x2"])
end

@testset "is_not_valid_variable_name" begin
    @syms x1::Bool T::Bool y::Bool
    @test !TransmissionMechanisms.is_not_valid_variable_name(x1)   # Returns false (valid variable name)
    @test TransmissionMechanisms.is_not_valid_variable_name(y)    # Returns true (invalid variable name)
    @test !TransmissionMechanisms.is_not_valid_variable_name(T)    # Returns false (valid variable name)
end

@testset "helper_sym_to_num" begin
    @syms x1::Bool x2::Bool x3::Bool T::Bool y::Bool
    nums = TransmissionMechanisms.helper_sym_to_num(x2)          # Returns [2]
    @test all(nums .== [2])
    nums = TransmissionMechanisms.helper_sym_to_num(T)          # Returns [nothing]
    @test all(isnothing.(nums))
    nums = TransmissionMechanisms.helper_sym_to_num([x1, x3]) # Returns [1, 3, nothing]
    @test all(nums .== [1, 3])
    @test_throws "Variable names should" TransmissionMechanisms.helper_sym_to_num(y)            # Error: Variable names should start with 'x' followed by a number.
end

@testset "get_terms" begin
    @syms x2::Bool x3::Bool Q(::Bool)::Real
    s = "x2 & !x3"
    cond = make_condition(s)
    term = TransmissionMechanisms.get_terms(cond)  # [Q(x2), -Q(x2 & x3)]
    @test all([isequal(x, y) for (x, y) in zip(term, [-Q(x3 & x2), Q(x2)])])
end

@testset "helper_Q" begin
    s = "x2 & !x3"
    cond = make_condition(s)
    terms, multiplier, variables, variable_nums = TransmissionMechanisms.helper_Q(cond)
    @test all(multiplier .== [1, -1])
    @test all(variable_nums .== [[2], [2, 3]])
end

@testset "contains_nots" begin
    s = "x2"
    term = make_condition(s)
    @test !TransmissionMechanisms.contains_nots(term)

    s = "x2 & !x3"
    term = make_condition(s)
    @test TransmissionMechanisms.contains_nots(term)
end