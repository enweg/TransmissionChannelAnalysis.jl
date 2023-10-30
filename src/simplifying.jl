import Base.&
import Base.|
import Base.!
import Base.string
import Base.show

"""
    Q(s::String)
    Q(s::String, m::Number)
    Q(s::Vector{String})
    Q(s::Vector{String}, m::Vector{Number})
    Q(i::Int)

Represents a transmission condition. 

We denote with Q(b), where b is a Boolean statement, a transmission question. 

**Important**: Each string in `s:` should be a Boolean statement involving
variables `x`` followed by a number, i.e. `x1`, `x2`, etc. Additionally, each
Boolean statement should contain only AND (&) and NOT (!) statements. 

**Important**: Users should only use the `Q(i::Int)` constructor. All other
constructors are for internal use only. Misuse of the other constructors easily
leads to mistakes. 

## Fields

- `vars::Vector{String}`: Contains the variables. These will be Boolean
  statements containing only AND and NOT.
- `multiplier::Vector{Number}`: Multiplier in front of Q(b). 

## Arguments

- `s::Union{String, Vector{String}}`: String representation of transmission
  condition. 
- `m::Union{Number, Vector{Number}}`: Multipliers for transmission conditions. 
- `i::Int`: Variable number.

## Usage

```julia

# Defining all variables in one go
x = [Q("x\$i") for i = 1:10]
q = (x[1] | x[2]) & !x[3]

# Alternatively variables can be defined separaterly
x1 = Q("x1")
x2 = Q("x2")
x3 = Q("x3")
q = (x1 | x2) & !x3

# The following are also valid
q = Q("x1 & !x3", 1)
q = Q(["x1", "x2", "x1 & x2"], [1, 1, -1])

# The following is NOT valid but does not yet throw an error or warning
q = Q("x1 | x2")  # DO NOT DO THIS!
```
"""
struct Q
    vars::Vector{String}
    multiplier::Vector{Number}
end
# TODO: Throw error if string is not valid.
Q(s::String) = Q([s], [1.0])
Q(s::String, m::Number) = Q([s], [m])
Q(s::Vector{String}) = Q(s, ones(size(s)))
Q(i::Int) = Q("x$i")

"""
    collect_terms(q::Q)

Collect all terms Q(b) for which the Boolean statement b is the same and sums
their multiplier. The result is a transmission condition for which each term
only appears ones, but with multipliers possibly different from plus-minus one. 

## Arguments

- `q::Q`: A transmission condition. See also [`Q`](@ref)

## Returns

- Another transmission condition of type `Q`. 

## Example

```julia
q = Q(["x1", "x1"], [1, 1])  
collect_terms(q)
# output: Q("x1", 2)
q = Q(["x1", "", "x1"], [1, 1, -1])  
collect_terms(q)
# output: Q("", 1)
```
"""
function collect_terms(q::Q)
    terms = Dict()
    for (v, m) in zip(q.vars, q.multiplier)
        m += get(terms, v, 0)
        terms[v] = m
    end
    vars = string.(keys(terms))
    mult = [terms[v] for v in vars]
    non_zero_mult = findall(!=(0), mult)
    vars = vars[non_zero_mult]
    mult = mult[non_zero_mult]
    return Q(vars, mult)
end

"""
    string_and(s1::String, s2::String)

Combine two strings using "&".
"""
function string_and(s1::String, s2::String)
    s2 == "" && return s1
    s1 == "" && return s2
    s = join([s1, s2], " & ")
    xs = sort(unique([m.match for m in eachmatch(r"(!{0,1}x\d+)", s)]); rev = true)
    s = join(xs, " & ")
    return s
end

"""
    check_contradiction(var_and::Vector{Int}, var_not::Vector{Int})
    check_contradiction(var_and::Vector{Vector{Int}}, var_not::Vector{Vector{Int}})
    
Check whether there is a contradiction of the form x1 & !x1. 

## Arguments

- `var_and::Union{Vector{Int}, Vector{Vector{Int}}}`: AND variable numbers
  obtained from [`get_varnums_and_multiplier`](@ref). 
- `var_not::Union{Vector{Int}, Vector{Vector{Int}}}`: NOT variable numbers
  obtained from [`get_varnums_and_multiplier`](@ref)

## Returns

1. `Bool` indicating whether there are any contradictions. 
2. `Vector{Bool}` indicating which elements yielded a contradiction.

## Notes

- This is used in [`remove_contradiction`](@ref) to remove contradicting terms.
  This speeds up simplification of terms, since the total number of terms can
  often be reduced. 

"""
function check_contradiction(var_and::Vector{Int}, var_not::Vector{Int})
    contradictions = [va in var_not for va in var_and]
    return any(contradictions), contradictions
end
function check_contradiction(var_and::Vector{Vector{Int}}, var_not::Vector{Vector{Int}})
    contradictions = [check_contradiction(va, vn)[1] for (va, vn) in zip(var_and, var_not)]
    return any(contradictions), contradictions
end

REMOVE_CONTRADICTIONS::Bool = true
"""
    remove_contradictions(q::Q)

Remove contradicting terms. 

A terms is deemed contradicting if it includes some "xi & !xi". This would
result in the entire Boolean statement to be false, and thus in the effect of
this terms to be zero. 

## Arguments

- `q::Q`: A transmission condition. See also [`Q`](@ref) and
  [`make_condition`](@ref). 

## Returns

- If `TransmissionMechanisms.REMOVE_CONTRADICTIONS == false`, then `q` will
  simply be returned again. 
- If `TransmissionMechanisms.REMOVE_CONTRADICTIONS == false`, then 
    1. If all terms are contradicting, the `Q("", 0)` will be retuned, which has
       a transmission effect of zero. 
    2. If some terms are non-contradicting, then a transmission condition
       consisting of only the non-contradicting terms will be returned. 

## Examples

```julia
TransmissionMechanisms.REMOVE_CONTRADICTIONS = true
q = TransmissionMechanisms.Q("x1", 1)
remove_contradictions(q)  # will return q again since no contradictions exist

q = TransmissionMechanisms.Q("x1 & !x1", 1)
remove_contradictions(q)  # Will return Q("", 0)

q = TransmissionMechanisms.Q(["x1 & !x1", "x1 & x2"], [1, 1])
remove_contradictions(q)  # Will return Q("x1 & x2", 1)
```
"""
function remove_contradictions(q::Q)
    global REMOVE_CONTRADICTIONS
    !REMOVE_CONTRADICTIONS && return q

    var_and, var_not, _ = get_varnums_and_multiplier(q)
    has_contradiction, contradictions =  check_contradiction(var_and, var_not)
    if has_contradiction
        all(contradictions) && return Q("", 0)  # This will later result in a zero
        return Q(q.vars[(!).(contradictions)], q.multiplier[(!).(contradictions)])
    end
    return q
end


"""
    Base.&(q1::Q, q2::Q)

Combine two transmission conditions using AND. 
"""
function (&)(q1::Q, q2::Q)
    if length(q1.vars) == 1
        q = collect_terms(Q([string_and(q2.vars[i], q1.vars[1]) for i = 1:lastindex(q2.vars)], q1.multiplier[1]*q2.multiplier))
        return remove_contradictions(q)
    elseif length(q2.vars) == 1
        q = collect_terms(Q([string_and(q1.vars[i], q2.vars[1]) for i = 1:lastindex(q1.vars)], q2.multiplier[1]*q1.multiplier))
        return remove_contradictions(q)
    else
        qs = [Q(q1.vars[i], q1.multiplier[i]) & q2 for i = 1:lastindex(q1.vars)]
        vars = reduce(vcat, [q.vars for q in qs])
        mults = reduce(vcat, [q.multiplier for q in qs])
        return collect_terms(Q(vars, mults))
    end
end

"""
    Base.|(q1::Q, q2::Q)

Combine two transmission conditions using OR. 
"""
function (|)(q1::Q, q2::Q)
    vars = q1.vars
    vars = vcat(vars, q2.vars)
    q = q1 & q2
    vars = vcat(vars, q.vars)
    mults = vcat(q1.multiplier, q2.multiplier)
    mults = vcat(mults, -1 * q.multiplier)
    return collect_terms(Q(vars, mults))
end

"""
    Base.!(q1::Q)

Return NOT the transmission condition if the condition involves more than one
variables. If the condition only involves one variables, then "!x1" is returned
where "1" is replaced by the respective variable number. 

Note: The decision not to simplify terms of the form "!x1" was made because
calculations usign the second calculation method in $WEGNER are
faster than having to simplify "!x1" type of terms and using the first
calculation method. 
"""
function (!)(q1::Q)
    if length(q1.vars) == 1 && count(r"x\d+", q1.vars[1]) == 1
        vars = "!$(q1.vars[1])"
        return Q(vars, q1.multiplier[1])
    end
    vars = vcat([""], q1.vars)
    mults = vcat([1.0], -1 * q1.multiplier)
    return collect_terms(Q(vars, mults))
end

"""
    Base.string(q::Q)

Obtain a string representation of a transmission condition.
"""
function Base.string(q::Q)
    s = []
    for (m, v) in zip(q.multiplier, q.vars)
        ms = string(m)
        if m % 1 == 0
            ms = string.(floor(Int, m))
        end
        if m == -1
            ms = "-"
        end
        if m == 1
            ms = ""
        end
        
        vs = string(v)
        if vs == ""
            vs = "T"
        end
            
        push!(s, "$(ms)Q($(vs))")
    end
    return s
end
function  Base.show(io::IO, ::MIME"text/plain", q::Q)
    if haskey(ENV, "SHOW_Q_AS_SUM") && ENV["SHOW_Q_AS_SUM"] == "true"
        s = join(string(q), " + ")
        return write(io, s)
    end
    s = join(string(q), "\n")
    return write(io, s)
end