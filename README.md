# TransmissionMechanisms.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://enweg.github.io/TransmissionMechanisms.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://enweg.github.io/TransmissionMechanisms.jl/dev/)
[![Build Status](https://github.com/enweg/TransmissionMechanisms.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/enweg/TransmissionMechanisms.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/enweg/TransmissionMechanisms.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/enweg/TransmissionMechanisms.jl)


Calculate the transmission effect along a collection of paths in an SVAR. See
paper for description. 

## Information for Users

Variables follow the convention that they start with "x" followed by a
number. The number corresponds to the variable number in the rewritten
structural form in the paper. For example, for a SVAR with $k$ variables, the
following mapping applies 

- $y_{1, t} \to x_1$
- $y_{2, t} \to x_2$
- ...
- $y_{k, t} \to x_k$
- $y_{1, t+1} \to x_{k+1}$ 
- ...

Transmission mechanisms are defined in terms of Boolean statements. To create
a transmission condition, first create a string corresponding to the desired
Boolean statement. For example, if the interest is in the transmission along
the paths that go through $x_2$ but not through $x_3$ or $x_4$,
then the string is of the form 

```julia
s = "x1 & !x3 & !x4".
```
This could equivalently be represented as `x1 & !(x3 | x4)`, however, the above
representation is often more efficient.  

String representations of Boolean statements can then be transformed into
transmission conditions/queries using `make_condition` as follows

```julia
cond = make_condition(s)
```

`cond` is now of type `Q`, which represents a transmission condition or
**Q**ery. The effect of this transmission query can be calculated using
`transmission` as follows

```julia
from_shock = 1
effect = transmission(from_shock, B, Qbb, cond)
effect = transmission(from_shock, B, Qbb, cond; method = :BQbb)  # equivalent to above.
```

where `B` and `Qbb` correspond to the structural matrices in the rewritten
structural form (also see the documentation for transmission). To obtain `B` and
`Qbb` from an estimated and (partially) identified SVAR, use
`to_structural_transmission_model`. If only IRFs are available, or if the IRF
method is preferred, then effect can be calculated by using `transmission` and
setting the keyword argument `method = :irfs` as follows

```julia
from_shock = 1
effect_irfs = transmission(from_shock, irfs, irfs_ortho, cond; method = :irfs)
```

where the `irfs` and `irfs_ortho` are the structural IRFs (with possibly `NaN`
columns indicating unidentified shocks) and `irfs_ortho` are the orthogonalised
IRFs. Both are in transmission form obtained using `to_transmission_irfs`. 

The `:BQbb` method is often the more efficient one, and as implemented here can
be used for all Boolean statements. Both the `:BQbb` and the `:irfs` method
return a vector with index `i` being the transmission effect on $x_i$. If $x_k$
is the variable with the highest subscript involved in the Boolean statement,
then the first `k` entries in the returned vector are `NaN` since interpretation
of those effects is nonsensical. 


## Internals

The internals of `TransmissionMechanisms.jl` all revolve around the type `Q`
which represents a transmission condition or query. It has two fields. The `vars`
field is a `Vector{String}` field with each element corresponding to a term in a
transmission condition. The second field is the `multiplier` field which is a
`Vector{Number}` and holds the multiplier for each term. For example

- The condition $Q[x_1]$ results in a single element in `vars` corresponding to
  `"x1"` and a single element in `multiplier` equal to `1`. 
- The condition $Q[x_1 \land x_2]$ results in a single element in `vars` corresponding
  to `"x2 & x1"` and a single element in `multiplier` corresponding to `1`. 
- The condition $Q[x_1 \lor x_2] = Q[x_1] + Q[x_2] - Q[x_1 \land x_2]$ results in three elements in `vars` corresponding
  to `"x1"`, `"x2"`, and `"x2 & x1"` and three elements in `multiplier`
  corresponding to `1`, `1`, and `-1` respectively.

Starting with simple queries of the form `Q("x1")`, and `Q("x2")`, the
overloaded operators yield simplifications consistent with the algebra in the
paper. This consistency cannot be guaranteed if one is not starting with these
"atomic" queries. 

`Q` follows three rules: 

1. Assume $Q[b] = \sum_{i = 1}^{N_1} m_i Q[b_i]$ and $Q[b^*] =\sum_{j=1}^{N_2}m^*_jQ[b^*_j]$. 
   Then $Q[b \land b^*] = \sum_{i=1}^{N_1}\sum_{j=1}^{N_2}m_im^*_jQ[b_i \land b^*_j]$
2. $Q[b \lor b^*] = Q[b] + Q[b^*] - Q[b \land b^*]$ 
3. $Q[b \land \neg b^*] = Q[b] - Q[b \land b^*]$

Rule (2) and (3) are included in the rules in the paper. Rule (1), on the other
hand, requires some more explanation. For this, I will first show that if $Q[b]
= \sum_{i = 1}^{N_1} m_i Q[b_i]$ then $Q[b \land b^*] = \sum_{i = 1}^{N_1} m_i
Q[b_i \land b^*]$.


In order for $Q[b]$ to simplify to $\sum_{i =1}^{N_1}Q[b_i]$, either rule (2) or
rule (3) must be applied at some point, because if $b$ was only to consist of
conjunctions, then the term would never be split (follows from the algebra).
Thus, if $b$ was only consisting of conjunctions, then the statement above is
trivially true. Now assume that simplification of $Q[b]$ to the sum of terms
takes $I$ iterations of applying rules (2) or (3) to some terms. At iteration
zero, where $Q[b]$ has not yet been split into other terms, $Q[b \land b^*] =
Q[b \land b^*]$ is trivially true. Thus, assume that the statement is true at
the beginning of some iteration $k$. We then must show that it is true at the
end of iteration $k$ and thus at the beginning of iteration $k+1$. The only two
operations that would further split a term are rules (2) and (3). Note though,
that 

1. $Q[(b_j \lor b_k) \land b^*] = Q[(b_j \land b^*) \lor (b_k \land b^*)] = Q[b_j \land b^*] + Q[b_k \land b^*] - Q[b_j \land b_k \land b^*]$
2. $Q[(b_j \land \neg b_k) \land b^*] = Q[b_j \land b^*] - Q[b_j \land b_k \land
   b^*]$

Thus, no-matter which rule is applied, the statement above would still hold for
the split term, and thus for the entire sum. This proofs the correctness by
induction and shows that $Q[b \land b^*] = \sum_{i=1}^{N_1}m_iQ[b_i \land b^*]$. 

To obtain the full statement, use the same arguments as above, but this time,
$b^*$ takes the role that $b$ took above, and $b_i$ takes the role that $b^*$
took above. This results in $Q[b_i \land b^*] = \sum_{j=1}^{N_2}m^*_jQ[b_i \land
b_j^*]$. Thus, $Q[b \land b^*] = \sum_{i=1}^{N_1}\sum_{j=1}^{N_2}m_im^*_jQ[b_i
\land b^*_j]$

We can thus implement the algebra in the paper by overwriting the AND (`&`), the
OR (`|`), and the NOT (`!`) operator for the type of `Q`, and implement the
rules above. This is what we internally do with the single exception that we do
not further simplify negations of single variables. Reason for this is
that it is faster to calculate the effect including the negation than it is to
simplify the negated statements to a state in which each term is only a
conjunction of ANDs. 

> There are obviously some more internal details. Please ask if you are
> interested. 