@syms Q(::Bool)::Real T::Bool

# Collect all terms in a transmission expression. For example, simplifying an
# expression might result in Q(x2 & !x3) = Q(x2) - Q(x2 & x3). Q(x2) and 
# -Q(x2 & x3) are collected.
collect_Qs_1 = @rule ((+)(~~y)) => ~~y
collect_Qs_2 = @rule (~x) => ~x
collect_Qs = Chain([collect_Qs_1, collect_Qs_2])

# Collect the multiplier infront of a term. 
collect_multiplier = @rule (~m * Q(~~y)) => ~m

# Q(a | b) = Q(a) + Q(b) - Q(a & b)
r_or_1 = @acrule (Q(~a | ~b)) => (Q(~a) + Q(~b) - Q(~a & ~b))
r_or_2 = @acrule (Q(~a | ~b) +(~~y)) => (Q(~a) + Q(~b) - Q(~a & ~b) + sum(~~y))
r_or_3 = @acrule (~m * Q(~a | ~b)) => (~m * Q(~a) + ~m * Q(~b) - ~m * Q(~a & ~b))
r_or_4 = @acrule (~m * Q(~a | ~b) +(~~y)) => (~m * Q(~a) + ~m * Q(~b) - ~m * Q(~a & ~b) + sum(~~y))
r_or = Chain([r_or_1, r_or_2, r_or_3, r_or_4])

# Q(a & !b) = Q(a) - Q(a & b)
r_and_not_1 = @acrule (Q(~a & !~b)) => (Q(~a) - Q(~a & ~b))
r_and_not_2 = @acrule (Q(~a & !~b) +(~~y)) => (Q(~a) - Q(~a & ~b) + sum(~~y))
r_and_not_3 = @acrule (~m * Q(~a & !~b)) => (~m * Q(~a) - ~m * Q(~a & ~b))
r_and_not_4 = @acrule (~m * Q(~a & !~b) +(~~y)) => (~m * Q(~a) - ~m * Q(~a & ~b) + sum(~~y))
r_and_not_5 = @acrule (Q(!~b & ~a)) => (Q(~a) - Q(~a & ~b))
r_and_not_6 = @acrule (Q(!~b & ~a) +(~~y)) => (Q(~a) - Q(~a & ~b) + sum(~~y))
r_and_not_7 = @acrule (~m * Q(!~b & ~a)) => (~m * Q(~a) - ~m * Q(~a & ~b))
r_and_not_8 = @acrule (~m * Q(!~b & ~a) +(~~y)) => (~m * Q(~a) - ~m * Q(~a & ~b) + sum(~~y))
r_and_not = Chain([
        r_and_not_1, r_and_not_2, r_and_not_3, r_and_not_4, 
        r_and_not_5, r_and_not_6, r_and_not_7, r_and_not_8
])

# Q(x & (y | z)) = Q((x & y) | (x & z))
r_or_cummutative_1 = @acrule (Q(~x & (~y | ~z))) => (Q((~x & ~y) | (~x & ~z)))
r_or_cummutative_2 = @acrule (Q(~x & (~y | ~z)) +(~~s)) => (Q((~x & ~y) | (~x & ~z)) + sum(~~s))
r_or_cummutative_3 = @acrule (~m * Q(~x & (~y | ~z))) => (~m * Q((~x & ~y) | (~x & ~z)))
r_or_cummutative_4 = @acrule (~m * Q(~x & (~y | ~z)) +(~~s)) => (~m * Q((~x & ~y) | (~x & ~z)) + sum(~~s))
r_or_cummutative_5 = @acrule (Q((~y | ~z) & ~x)) => (Q((~x & ~y) | (~x & ~z)))
r_or_cummutative_6 = @acrule (Q((~y | ~z) & ~x &(~~s))) => (Q((~x & ~y &(~~s...)) | (~x & ~z &(~~s...))))
r_or_cummutative_7 = @acrule (~m * Q((~y | ~z) & ~x)) => (~m * Q((~x & ~y) | (~x & ~z)))
r_or_cummutative_8 = @acrule (~m * Q((~y | ~z) & ~x &(~~s))) => (~m * Q((~x & ~y &(~~s...)) | (~x & ~z &(~~s...))))
r_or_cummutative = Chain([
        r_or_cummutative_1, r_or_cummutative_2, r_or_cummutative_3, 
        r_or_cummutative_4, r_or_cummutative_5, r_or_cummutative_6,
        r_or_cummutative_7, r_or_cummutative_8
])

# Q(x & (x | y)) = Q(x)
r_absorption_and_1 = @acrule (Q(~x & (~x | ~y))) => (Q(~x))
r_absorption_and_2 = @acrule (Q(~x & (~x | ~y)) +(~~y)) => (Q(~x) + sum(~~y))
r_absorption_and_3 = @acrule (~m * Q(~x & (~x | ~y))) => (~m * Q(~x))
r_absorption_and_4 = @acrule (~m * Q(~x & (~x | ~y)) +(~~y)) => (~m * Q(~x) + sum(~~y))
r_absorption_and = Chain([r_absorption_and_1, r_absorption_and_2, r_absorption_and_3, r_absorption_and_4])

# Q(x | (x & y)) = Q(x)
r_absorption_or_1 = @acrule (Q(~x | (~x & ~y))) => (Q(~x))
r_absorption_or_2 = @acrule (Q(~x | (~x & ~y)) +(~~s)) => (Q(~x) + sum(~~s))
r_absorption_or_3 = @acrule (~m * Q(~x | (~x & ~y))) => (~m * Q(~x))
r_absorption_or_4 = @acrule (~m * Q(~x | (~x & ~y)) +(~~s)) => (~m * Q(~x) + sum(~~s))
r_absorption_or = Chain([r_absorption_or_1, r_absorption_or_2, r_absorption_or_3, r_absorption_or_4])

# Q(x | (y & z)) = Q((x | y) & (x | z))
r_and_cummutative_1 = @acrule (Q(~x | (~y & ~z))) => (Q((~x | ~y) & (~x | ~z)))
r_and_cummutative_2 = @acrule (Q(~x | (~y & ~z)) +(~~s)) => (Q((~x | ~y) & (~x | ~z)) + sum(~~s))
r_and_cummutative_3 = @acrule (~m * Q(~x | (~y & ~z))) => (~m * Q((~x | ~y) & (~x | ~z)))
r_and_cummutative_4 = @acrule (~m * Q(~x | (~y & ~z)) +(~~s)) => (~m * Q((~x | ~y) & (~x | ~z)) + sum(~~s))
r_and_cummutative = Chain([r_and_cummutative_1, r_and_cummutative_2, r_and_cummutative_3, r_and_cummutative_4])

# Q(!(x | y)) = Q(!x & !y)
r_not_or_1 = @acrule (Q(!(~x | ~y))) => (Q(!(~x) & !(~y)))
r_not_or_2 = @acrule (Q(!(~x | ~y)) +(~~s)) => (Q(!(~x) & !(~y)) + sum(~~s))
r_not_or_3 = @acrule (~m * Q(!(~x | ~y))) => (~m * Q(!(~x) & !(~y)))
r_not_or_4 = @acrule (~m * Q(!(~x | ~y)) +(~~s)) => (~m * Q(!(~x) & !(~y)) + sum(~~s))
r_not_or = Chain([r_not_or_1, r_not_or_2, r_not_or_3, r_not_or_4])

# Q(!(x & y)) = Q(!x | !y)
r_not_and_1 = @acrule (Q(!(~x & ~y))) => (Q(!(~x) | !(~y)))
r_not_and_2 = @acrule (Q(!(~x & ~y)) +(~~s)) => (Q(!(~x) | !(~y)) + sum(~~s))
r_not_and_3 = @acrule (~m * Q(!(~x & ~y))) => (~m * Q(!(~x) | !(~y)))
r_not_and_4 = @acrule (~m * Q(!(~x & ~y)) +(~~s)) => (~m * Q(!(~x) | !(~y)) + sum(~~s))
r_not_and = Chain([r_not_and_1, r_not_and_2, r_not_and_3, r_not_and_4])

# Q(!x) = Q - Q(x)
r_not_1 = @acrule (Q(!~x)) => (Q(T) - Q(~x))
r_not_2 = @acrule (Q(!~x) +(~~s)) => (Q(T) - Q(~x) + sum(~~s))
r_not_3 = @acrule (~m * Q(!~x)) => (~m * Q(T) - ~m * Q(~x))
r_not_4 = @acrule (~m * Q(!~x) +(~~s)) => (~m * Q(T) - ~m * Q(~x) + sum(~~s))
r_not = Chain([r_not_1, r_not_2, r_not_3, r_not_4])

# reshuffle ANDs such that a not is in front. This helps the simplification process. 
r_reshuffle_not_1 = @acrule (Q((~x & !~y) &(~~z))) => (Q((&)(vcat(~x, ~~z)...) & !(~y)))
r_reshuffle_not_2 = @acrule (Q((~x & !~y) &(~~z)) +(~~s)) => (Q((&)(vcat(~x, ~~z)...) & !(~y)) + sum(~~s))
r_reshuffle_not_3 = @acrule (~m * Q((~x & !~y) &(~~z))) => (~m * Q((&)(vcat(~x, ~~z)...) & !(~y)))
r_reshuffle_not_4 = @acrule (~m * Q((~x & !~y) &(~~z)) +(~~s)) => (~m * Q((&)(vcat(~x, ~~z)...) & !(~y)) + sum(~~s))
r_reshuffle_not = Chain([r_reshuffle_not_1, r_reshuffle_not_2, r_reshuffle_not_3, r_reshuffle_not_4])

# Collect all the `&` terms. 
r_collect_ands_1 = @acrule (Q((&)(~~y))) => (Q((&)(collect_and_terms(~~y)...)))
r_collect_ands_2 = @acrule (Q((&)(~~y) +(~~s))) => (Q((&)(collect_and_terms(~~y)...) + sum(~~s)))
r_collect_ands_3 = @acrule (~m * Q((&)(~~y))) => (~m * Q((&)(collect_and_terms(~~y)...)))
r_collect_ands_4 = @acrule (~m * Q((&)(~~y) +(~~s))) => (~m * Q((&)(collect_and_terms(~~y)...) + sum(~~s)))
r_collect_ands = Chain([r_collect_ands_1, r_collect_ands_2, r_collect_ands_3, r_collect_ands_4])

ch = Chain([
            r_or, 
            r_and_not, 
            r_or_cummutative, 
            r_absorption_and, 
            r_absorption_or, 
            r_not_or, 
            r_not_and, 
            r_not, 
            r_reshuffle_not, 
            r_collect_ands
]);

# Helper rules to collect all unique AND terms. 
r_ands_1 = @acrule (Q(~x &(~~y))) => (x = ~x, y = ~~y)
r_ands_2 = @acrule (~x &(~~y)) => (x = ~x, y = ~~y)
r_ands_3 = @acrule (~m * Q(~x &(~~y))) => (x = ~x, y = ~~y)
r_ands_4 = @acrule (Q(~x)) => (x = ~x, y = [~x])
r_ands_5 = @acrule (~m * Q(~x)) => (x = ~x, y = [~x])
r_ands = Chain([r_ands_1, r_ands_2, r_ands_3, r_ands_4, r_ands_5])