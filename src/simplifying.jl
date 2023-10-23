"""
    collect_and_terms(expr::Union{SymbolicUtils.BasicSymbolic{Bool}, Vector{SymbolicUtils.BasicSymbolic{Bool}}}; rev=true)

Collect unique AND terms from a boolean expression or a vector of boolean expressions.

This function extracts unique AND terms from the input boolean expression(s) and 
returns them as a vector of `SymbolicUtils.BasicSymbolic{Bool}`. The input `expr` 
can be either a single boolean expression or a vector of boolean expressions.

## Arguments

- `expr::Union{SymbolicUtils.BasicSymbolic{Bool}, Vector{SymbolicUtils.BasicSymbolic{Bool}}}`: 
  A single boolean expression or a vector of boolean expressions from which AND 
  terms will be collected.
- `rev::Bool=true`: (Optional) If `true`, the resulting AND terms are sorted in 
  descending order based on their string representation. If `false`, the sorting 
  is in ascending order.

## Returns

- Returns a vector of `SymbolicUtils.BasicSymbolic{Bool}` containing unique AND terms 
  extracted from the input expression(s).

## Examples

```julia
s = "x2 & x3 & !x2 & x5"
cond = make_condition(s)
collect_and_terms(cond)
````

"""
function collect_and_terms(expr::SymbolicUtils.BasicSymbolic{Bool}; rev = true)
    stack = Any[expr]
    ands = []
    while length(stack) > 0
        ex = pop!(stack)
        out = r_ands(ex)
        if isequal(out, ex)
            push!(ands, out)
        else
            push!(stack, out.x)
            for y in out.y
                push!(stack, y)
            end
        end
    end
    ands_str = string.(ands)
    ps = sortperm(ands_str; rev = rev)
    return unique(ands[ps])
end
function collect_and_terms(expr::Vector{SymbolicUtils.BasicSymbolic{Bool}}; rev = true)
    out = [collect_and_terms(ex) for ex in expr]
    out = vcat(out...)
    out_str = string.(out)
    ps = sortperm(out_str; rev = rev)
    return out[ps]
end


"""
    is_not_valid_variable_name(sym::SymbolicUtils.BasicSymbolic{Bool})

Check if a given symbol is a valid variable name in the context of transmission
mechanisms.

This function validates whether the provided symbol `sym` is a valid variable
name for use in transmission mechanisms and in
[`create_transmission_function`](@ref). In the context of transmission
mechanisms, valid variable names are of the form 'x' followed by one or more
digits, or the symbol 'T' representing `true`.

## Arguments

- `sym::SymbolicUtils.BasicSymbolic{Bool}`: Symbol to be checked for validity as
  a variable name in transmission mechanisms.

## Returns

- Returns `true` if the input symbol is **not** a valid variable name, and `false`
  otherwise.

## Examples

```julia
@syms x1::Bool T::Bool y::Bool
is_not_valid_variable_name(:x1)   # Returns false (valid variable name)
is_not_valid_variable_name(:y)    # Returns true (invalid variable name)
is_not_valid_variable_name(:T)    # Returns false (valid variable name)

"""
function is_not_valid_variable_name(sym::SymbolicUtils.BasicSymbolic{Bool})
    sym_str = string(sym)
    if !isnothing(match(r"^x\d+$|^T$", sym_str))
        return false
    end
    return true
end


"""
    helper_sym_to_num(sym_vec::Union{SymbolicUtils.BasicSymbolic{Bool}, Vector{SymbolicUtils.BasicSymbolic{Bool}}})

Convert symbolic variables to numerical indices.

This function takes a single boolean symbolic variable or an array of boolean
symbolic variables `sym_vec` and converts them into numerical indices. Each
index corresponds to the variable index in the structural transmission model and
thus to an index in the IRF matrix. Valid variable names should start with 'x'
followed by a number or be the symbol 'T'. If invalid variable names are
encountered, an error is raised.

## Arguments

- `sym_vec::Union{SymbolicUtils.BasicSymbolic{Bool}, Vector{SymbolicUtils.BasicSymbolic{Bool}}}`: 
  A single boolean symbolic variable or an array of boolean symbolic variables 
  representing variables to be converted to numerical indices.

## Returns

- Returns a sorted vector of numerical indices corresponding to the input symbols. 
  If a symbol is 'T', it is represented as `nothing` in the output vector.

## Examples

```julia
@syms x1::Bool x2::Bool x3::Bool T::Bool y::Bool
helper_sym_to_num(:x2)          # Returns [2]
helper_sym_to_num([:x1, :x3, :T]) # Returns [1, 3, nothing]
helper_sym_to_num(:y)            # Error: Variable names should start with 'x' followed by a number.

"""
function helper_sym_to_num(sym_vec::Union{SymbolicUtils.BasicSymbolic{Bool}, Vector{SymbolicUtils.BasicSymbolic{Bool}}})
    if !isa(sym_vec, AbstractVector)
        sym_vec = [sym_vec]
    end
    any(is_not_valid_variable_name.(sym_vec)) && error("Variable names should start with 'x' followed by a number. See the paper for details.")
    rx = r"[0-9]*$|T$"
    out = [match(rx, string(sym)).match for sym in sym_vec]
    out = [o == "T" ? nothing : parse(Int, o) for o in out]
    return sort(out)
end

"""
    get_terms(condition::SymbolicUtils.BasicSymbolic{Bool})

Simplify a Boolean condition into a sum of only conjunctions. 

Given a valid Boolean condition in the context of transmission mechanisms, the
condition is recursively simplified until the result only consists of additive
terms with each remaining condition being purely a conjunction. For example, The
condition `x2 & !x3` corresponds to the transmission statement `Q(x2 & !x3)`
which is simplified to `Q(x2) - Q(x2 & x3)`. The function will then return
`[Q(x2), -Q(x2 & x3)]`. Each term can then be easily evaluated using IRFs.  

## Arguments

- `condition::SymbolicUtils.BasicSymbolic{Bool}`: A symbolic boolean condition.
  Can be created using [`make_condition`](@ref). Valid variable names are `x`
  followed by a number, or `T` respresenting true. All transmission paths must
  satisfy this boolean statement. 

## Returns 

- Returns a vector of terms, each of which consists of a conjunction. 

## Examples

```julia
s = "x2 & !x3"
cond = make_condition(s)
term = get_terms(cond)  # [Q(x2), -Q(x2 & x3)]
```

"""
function get_terms(condition::SymbolicUtils.BasicSymbolic{Bool})
    ex = Q(condition)
    stack = Any[ex]
    terms = Any[]
    while length(stack) > 0
        ex = pop!(stack)
        simplified = ch(ex)
        if isequal(ex, simplified)
            push!(terms, ex)
        else
            simplified_terms = collect_Qs(simplified)
            if isnothing(simplified_terms)
                push!(stack, simplified)
                continue
            end
            if !isa(simplified_terms, AbstractVector)
                simplified_terms = [simplified_terms]
            end
            for st in simplified_terms
                push!(stack, st)
            end
        end
    end
    return terms
end


"""
    helper_Q(condition::SymbolicUtils.BasicSymbolic{Bool})

Simplifies the valid Boolean `condition` in the context of transmission
mechanisms and returns all necessary information to evaluate the transmission
effect using information in Impulse Response Functions (IRFs). 

## Arguments 

- `condition::SymbolicUtils.BasicSymbolic{Bool}`: A valid symbolic Boolean
  condition. Can be created using [`make_condition`](@ref). 
  
## Returns 

1. A vector containing all the additive terms in the simplified representation.
   Each term only involves conjunctions and can thus be calculated using
   information in IRFs. 
2. A vector of multipliers applied to each term. 
3. A vector of variables involved in the transmission query.
4. A vector of vectors. Each inner vector corresponds to a term and provides the
   variable numbers involved in the term. These, together with the multipliers,
   can be used to calculate the transmission effect. 

"""
function helper_Q(condition::SymbolicUtils.BasicSymbolic{Bool})
    terms = get_terms(condition)
    term = (+)(terms...)
    terms = collect_Qs(term)
    multiplier = collect_multiplier.(terms)
    if !isa(multiplier, AbstractVector)
        multiplier = [multiplier]
    end
    multiplier = [isnothing(m) ? 1 : m for m in multiplier]
    variables = collect_and_terms.(terms; rev = false)
    variable_nums = helper_sym_to_num.(variables)
    return terms, multiplier, variables, variable_nums
end


"""
    contains_nots(term::SymbolicUtils.BasicSymbolic{Bool})

Check whether a `term` still contains a NOT statement. 

To be able to calculate transmission effects using IRFs, Boolean conditions need
to be simplified to a state in which they only involve ANDs in all terms. Thus,
this function checks whether something went wrong in the simplification process. 

## Arguments

- `term::SymbolicUtils.BasicSymbolic{Bool}`: A `term` as returned by
  [`helper_Q`](@ref) or [`get_terms`](@ref). 

## Returns

- Returns `true` if the `term` includes any NOT and `false` otherwise. 

"""
function contains_nots(term::SymbolicUtils.BasicSymbolic{Bool})
    term_string = string(term)
    return contains(term_string, "!")
end