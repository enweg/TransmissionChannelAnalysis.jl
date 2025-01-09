@testset "map y to x and back" begin
    order = [3, 1, 2]
    s_y = "y_{1,2} & y_{3,1}"
    s_x = map_y_to_x(s_y, order)  # Returns: "x8 & x4"
    @test s_y == map_x_to_y(s_x, order)

    # introducing spaces
    s_y2 = "y_{1, 2} & y_{3, 1}"
    s_x = map_y_to_x(s_y2, order)  # Returns: "x8 & x4"
    # compare to s_y because function returns without spaces
    @test s_y == map_x_to_y(s_x, order)

    # introducing not
    order = [3, 1, 2]
    s_y = "!y_{1,2} & y_{3,1}"
    s_x = map_y_to_x(s_y, order)  # Returns: "x8 & x4"
    @test s_y == map_x_to_y(s_x, order)

    # introducing zeros
    order = [3, 1, 2]
    s_y = "!y_{1,0} & y_{3,1}"
    s_x = map_y_to_x(s_y, order)  # Returns: "x8 & x4"
    @test s_y == map_x_to_y(s_x, order)

    # introducing parentheses
    order = [3, 1, 2]
    s_y = "!(y_{1,0} & y_{3,1})"
    s_x = map_y_to_x(s_y, order)  # Returns: "x8 & x4"
    @test s_y == map_x_to_y(s_x, order)
end

@testset "slide_in!" begin
    A = reduce(hcat, fill(i, 3, 3) for i = 3:-1:1)
    B = zeros(12, 12)
    TransmissionChannelAnalysis.slide_in!(B, A)
    @test all(B[1:3, 1:3] .== 1)
    @test all(B[1:3, 4:end] .== 0)
    @test all(B[4:6, 1:3] .== 2)
    @test all(B[4:6, 4:6] .== 1)
    @test all(B[4:6, 7:end] .== 0)
    @test all(B[7:9, 1:3] .== 3)
    @test all(B[7:9, 4:6] .== 2)
    @test all(B[7:9, 7:9] .== 1)
    @test all(B[7:9, 10:end] .== 0)
    @test all(B[10:end, 1:3] .== 0)
    @test all(B[10:end, 4:6] .== 3)
    @test all(B[10:end, 7:9] .== 2)
    @test all(B[10:end, 10:end] .== 1)

    @test_throws "A cannot slide" TransmissionChannelAnalysis.slide_in!(A, zeros(11, 12))
end
