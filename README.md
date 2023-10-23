# TransmissionMechanisms.jl

Calculate the transmission effect along a collection of paths in an SVAR. See
paper for description. 

## Example

Say you have a SVAR(4) with five endogenous variables. You are interested in the
effect going through $y_{3}$ in any period. Since there are five endogenous
variables, we have the following variable definitions: 

1. $y_{3, 0} \to x_3$
2. $y_{3, 1} \to x_8$
3. $y_{3, h} \to x_{3+5h}$

We are therefore interested in the transmission effect described by the Boolean
statement $\lor_{i=0}^h x_{3 + 5i}$ which can be calculated in the following
way: 

```julia
using TransmissionMechanisms

from = 1
h = 3
s = join(["x$(3 + 5*i)" for i = 0:h], " | ")
cond = make_condition(s)
tf = create_transmission_function(1, cond)
tf(irfs, irfs_ortho)
```

In the code above, we assumed that the structural shock of interest is the
first, and we assumed that `irfs` and `irfs_ortho` have already been created.
For details see the documentation. 

> Any Boolean statement can be used, as long as the variables start with $x$
> followed by a number. This follows the paper. AND is represented by `&`, NOT
> is represented by `!`, and OR is represented by `|`. 