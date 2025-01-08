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
