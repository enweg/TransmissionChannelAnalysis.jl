@doc raw"""

    make_structural_B(A0_ortho::AbstractMatrix, Aplus_ortho::AbstractMatrix, p::Int, max_horizon::Int)

Create the `B` matrix in the structural representation in $WEGNER. 

Given a SVAR(p) in the form of 
```math
y_t'A_0 = y_{t-1}'A_1 + \dots + y_{t-p}'A_p + \varepsilon_t' = [y_{t-1}', \dots, y_{t-p}']A^+ + \varepsilon_t',
```
and a pre-specified maximum horizon of ``h``, the SVAR(p) over the next `h`
horizons can be written as 
```math
x = Bx + \mathbb{Q}\epsilon,
```
with the structure of ``B`` and ``\mathbb{Q}`` given in the paper. 

## Arguments

- `A0_ortho::AbstractMatrix`: The contemporaneous matrix of the SVAR(p) obtained
  using a Cholesky decomposition. **Note that this would be an upper-triangular
  matrix in the representation chosen here.** 
- `Aplus_ortho::AbstractMatrix`: The structural lag matrix of the SVAR(p)
  obtained using a Cholesky decomposition. 
- `p::Int`: The order of the SVAR(p). 
- `max_horizon::Int`: The maximum horizon to consider for the transmission
  model. 

## Returns

- Returns ``B`` in the above structural representation of the SVAR(p) over the
  next ``h`` horizons. 

## Notes

- See also [`make_structural_Omega`](@ref) and [`to_structural_transmission_model`](@ref).

"""
function make_structural_B(A0_ortho::AbstractMatrix, Aplus_ortho::AbstractMatrix, p::Int, max_horizon::Int)
    _, k = size(A0_ortho)

    L = A0_ortho'
    D = inv(diagm(diag(L)))
    DL = D*L
    B = D*Aplus_ortho'
    
    B = reduce(hcat, [B[:, ((i-1)*k+1):(i*k)] for i=p:-1:1])
    B = hcat(B, I-DL)
    Bs = [hcat(
            zeros(k, max(0, (i+1)*k - k*(p+1))),
            B[:, max(1, (end-(i+1)*k+1)):end], 
            zeros(k, (max_horizon+1)*k - (i+1)*k)
          ) for i in 0:max_horizon]
    B = reduce(vcat, Bs)
    return B
end

@doc raw"""

    make_structural_B(A0_ortho::AbstractMatrix, Aplus_ortho::AbstractMatrix, p::Int, max_horizon::Int)

Create the `B` matrix in the structural representation in $WEGNER. 

Given a SVAR(p) in the form of 
```math
y_t'A_0 = y_{t-1}'A_1 + \dots + y_{t-p}'A_p + \varepsilon_t' = [y_{t-1}', \dots, y_{t-p}']A^+ + \varepsilon_t',
```
and a pre-specified maximum horizon of ``h``, the SVAR(p) over the next `h`
horizons can be written as 
```math
x = Bx + \mathbb{Q}\epsilon,
```
with the structure of ``B`` and ``\mathbb{Q}`` given in the paper. 

## Arguments

- `irf0::AbstractMatrix`: The identified impact responses. If the `i`th shock
  has not been identified, then the `i`th column consists of `NaN`. Matrix
  should be of dimension `k` times `k` with `k` being the number of variables in
  the SVAR(p).
- `A0_ortho::AbstractMatrix`: The contemporaneous matrix of the SVAR(p) obtained
  using a Cholesky decomposition. **Note that this would be an upper-triangular
  matrix in the representation chosen here.** 
- `max_horizon::Int`: The maximum horizon to consider for the transmission
  model. 

## Returns

- Returns ``\mathbb{Q}`` in the above structural representation of the SVAR(p) over the
  next ``h`` horizons. 

## Notes

- See also [`make_structural_B`](@ref) and [`to_structural_transmission_model`](@ref).

"""
function make_structural_Omega(irf0::AbstractMatrix, A0_ortho::AbstractMatrix, max_horizon::Int)
    
    m, k = size(A0_ortho)
    
    L = A0_ortho'
    D = inv(diagm(diag(L)))
    
    Qt = L*irf0
    DQt = D*Qt
    eye_h = diagm(ones(max_horizon+1))
    Omega = kron(eye_h, DQt)
    
    return Omega
end

@doc raw"""

    make_structural_B(A0_ortho::AbstractMatrix, Aplus_ortho::AbstractMatrix, p::Int, max_horizon::Int)

Create the `B` matrix in the structural representation in $WEGNER. 

Given a SVAR(p) in the form of 
```math
y_t'A_0 = y_{t-1}'A_1 + \dots + y_{t-p}'A_p + \varepsilon_t' = [y_{t-1}', \dots, y_{t-p}']A^+ + \varepsilon_t',
```
and a pre-specified maximum horizon of ``h``, the SVAR(p) over the next `h`
horizons can be written as 
```math
x = Bx + \mathbb{Q}\epsilon,
```
with the structure of ``B`` and ``\mathbb{Q}`` given in the paper. 

## Arguments

- `irf0::AbstractMatrix`: The identified impact responses. If the `i`th shock
  has not been identified, then the `i`th column consists of `NaN`. Matrix
  should be of dimension `k` times `k` with `k` being the number of variables in
  the SVAR(p).
- `A0_ortho::AbstractMatrix`: The contemporaneous matrix of the SVAR(p) obtained
  using a Cholesky decomposition. **Note that this would be an upper-triangular
  matrix in the representation chosen here.**
- `Aplus_ortho::AbstractMatrix`: The structural lag matrix of the SVAR(p)
  obtained using a Cholesky decomposition. 
- `p::Int`: The order of the SVAR(p). 
- `max_horizon::Int`: The maximum horizon to consider for the transmission
  model. 

## Returns

- Returns ``B`` and ``\mathbb{Q}`` in the above structural representation of the SVAR(p) over the
  next ``h`` horizons in this order. 

"""
function to_structural_transmission_model(irf0::AbstractMatrix, A0_ortho::AbstractMatrix, Aplus_ortho::AbstractMatrix, p::Int, max_horizon::Int)
    k = size(A0_ortho, 1)
    n_ex = size(Aplus_ortho, 1) - p*k
    B = make_structural_B(A0_ortho, Aplus_ortho[(n_ex+1):end, :], p, max_horizon)
    Omega = make_structural_Omega(irf0, A0_ortho, max_horizon)
    return B, Omega
end
