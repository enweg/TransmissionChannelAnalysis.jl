
"""
    map_y_to_x(i::Int, t::Int, K::Int, order::AbstractVector{<:Int}) -> Int

Maps a variable in the original transmission condition (using variables `y`) 
to its corresponding systems form (using variables `x`).

# Arguments
- `i::Int`: Index of the variable in the original system.
- `t::Int`: Time index of the variable in the original system.
- `K::Int`: Total number of variables in the system.
- `order::AbstractVector{<:Int}`: Vector defining the ordering of variables 
  (defined by the transmission matrix).

# Returns
- An integer representing the index of the corresponding `x` variable in the systems form.

# Notes
- The contemporaneous period is denoted is period `0`.
"""
map_y_to_x(i::Int, t::Int, K::Int, order::AbstractVector{<:Int}) = K * t + findfirst(==(i), order)

"""
    map_y_to_x(s_y::String, order::AbstractVector{<:Int}) -> String

Transforms a string representing a transmission condition in terms of the original
variables `y` to a string using the systems form variables `x` based on a specified
ordering.

# Arguments
- `s_y::String`: A string representing the transmission condition using original
  variables `y`, e.g. "y_{1,0} & !y_{2,0}".
- `order::AbstractVector{<:Int}`: Ordering defined in the transmission matrix.

# Returns
- A string representing the same transmission question but using variables 
  of the  systems form (`x`) rather than those of the dynamic (original) system (`y`).

# Example
```julia
order = [3, 1, 2]
s_y = "y_{1, 2} & y_{3, 1}"
map_y_to_x(s_y, order)  # Returns: "x8 & x4"
```
"""
function map_y_to_x(s_y::String, order::AbstractVector{<:Int})
    s_x = s_y
    for m in eachmatch(r"y_{(\d+),\s*(\d+)}", s_y)
        i, t = parse.(Int64, m.captures)
        xi = map_y_to_x(i, t, length(order), order)
        s_x = replace(s_x, m.match => "x$(xi)")
    end
    return s_x
end

"""
    map_x_to_y(xi::Int, order::AbstractVector{<:Int}) -> Tuple{Int, Int}

Maps a variable from the systems form (using variables `x`) back to its original
transmission condition form (using variables `y`). This function is the reverse of
`map_y_to_x`.

# Arguments
- `xi::Int`: Index of the variable in the systems form.
- `order::AbstractVector{<:Int}`: Ordering defined in the transmission matrix.

# Returns
- A tuple `(i, t)` where `i` is the index of the variable in the original system
  and `t` is its time index.

"""
function map_x_to_y(xi::Int, order::AbstractVector{<:Int})
    K = length(order)
    t = floor(Int, (xi - 1) / K)
    i = xi - t * K  # under the current transmission matrix
    return order[i], t  # under the original order
end

"""
    map_x_to_y(s_x::String, order::AbstractVector{<:Int}) -> String

Transforms a string representing a transmission condition in terms of the systems
form variables `x` back to a string using the original variables `y`. This function
is the reverse of `map_y_to_x`.

# Arguments
- `s_x::String`: A string representing the transmission condition using systems
  form variables `x`.
- `order::AbstractVector{<:Int}`: Ordering defined by the transmission matrix.

# Returns
- A string representing the same transmission question but using variables 
  of the dynamic (original) system (`y`) rather than those of the  systems form (`x`).

# Example
```julia
order = [3, 1, 2]
s_x = "x8 & x4"
map_x_to_y(s_x, order)  # Returns: "y_{1, 2} + y_{3, 1}"
```
"""
function map_x_to_y(s_x::String, order::AbstractVector{<:Int})
    s_y = s_x
    for m in eachmatch(r"x(\d+)", s_x)
        xi = parse(Int64, m.captures[1])
        i, t = map_x_to_y(xi, order)
        s_y = replace(s_y, m.match => "y_{$i,$t}")
    end
    return s_y
end

"""
Create a permutation matrix.
"""
permmatrix(order::AbstractVector{<:Int}) = I(length(order))[order, :]

"""
Slides A into B. 

# Notes
- used to construct the B and Omega matrices. 

"""
function slide_in!(B::AbstractMatrix, A::AbstractMatrix)
    mod(size(B, 1), size(A, 1)) == 0 || error("A cannot slide into B because the number of rows of B is not an integer multiple of the number of rows of A.")

    size(B, 1) == size(B, 2) || error("B must be square.")

    n_horizontal_blocks = floor(Int, size(B, 1) / size(A, 1))
    K = size(A, 1)
    for i = 1:n_horizontal_blocks
        block = A[:, max(1, (end - i * K + 1)):end]
        r = ((i-1)*K+1):(i*K)
        c = (i*K-size(block, 2)+1):(i*K)
        B[r, c] .= block
    end
    return B
end

"""
    set_remove_contradictions(b::Bool)

Sets REMOVE_CONTRADICTIONS to `b`. 

If `REMOVE_CONTRADICTIONS==true`, then contradicting terms are removed before
the transmission effect is computed. A term is deemed contradicting if it 
includes some "xi & !xi". This would result in the entire Boolean statement to 
be false, and thus in the effect of this terms to be zero. 
"""
function set_remove_contradictions(b::Bool)
    REMOVE_CONTRADICTIONS[] = b
end

"""
    ql(A::AbstractMatrix) -> Q, L

Compute the QL decomposition of a matrix `A`, such that `A = Q * L`, where:

- `Q` is an orthogonal matrix (`Q'Q = I`),
- `L` is a lower-triangular matrix with positive diagonal entries.

# Arguments
- `A::AbstractMatrix`: The input matrix to decompose.

# Returns
- `Q`: An orthogonal matrix of the same size as `A`, representing the "Q" in 
  the QL decomposition.
- `L`: A lower-triangular matrix with positive diagonal elements.
"""
function ql(A::AbstractMatrix)
    k = size(A, 2)
    P = permmatrix(collect(k:-1:1))
    qr_decomposition = qr(A * P)
    L = P * qr_decomposition.R * P
    Q = qr_decomposition.Q * P
    S = diagm(sign.(diag(L)))
    L = S*L
    Q = Q * S
    return Q, L
end


