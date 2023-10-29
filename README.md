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
the paths that go through $x_2$ but not through $x_3$ or $x_4$ onto $x_5$,
then the string is of the form 

```julia
s = "x1 & !x3 & !x4".
```
This could equivalently be represented as `x1 & !(x3 | x4)`, however, the above
representation is more efficient.  

String representations of Boolean statements can then be transformed into
transmission conditions using `make_condition` as follows

```julia
cond = make_condition(s)
```

`cond` is now of type `Q`, which represents a transmission condition or
**Q**ery. The effect of this transmission query can be calculated using
`transmission` as follows

```julia
from_shock = 1
effect = transmission(from_shock, B, Qbb, cond)
```

where `B` and `Qbb` correspond to the structural matrices in the rewritten
structural form (also see the documentation for transmission). To obtain `B` and
`Qbb` from an estimated and (partially) identified SVAR, use
`to_structural_transmission_model`. 
