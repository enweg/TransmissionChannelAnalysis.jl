# Simulate a structural system with `k` variables and a single horizon.
# Starting from a Cholesky factor `L` (lower triangular, positive diagonal) and
# a rotation `Q`, the structural IRFs are `L * Q'` and the orthogonalised IRFs
# are `L`. The systems form is then `B = I - D * inv(L)` (strictly lower
# triangular) and `Omega = D * Q'` with `D = Diagonal(L)`, so that
# `(I - B) \ Omega == L * Q'`. This allows to compare the :BOmega and the :irfs
# method on the same system.
function simulate_masking_system(k::Int; seed::Int=6150533)
    rng = Xoshiro(seed)
    Q = Matrix(qr(randn(rng, k, k)).Q)
    L = tril(randn(rng, k, k), -1) + Diagonal(1 .+ rand(rng, k))
    D = Diagonal(diag(L))
    B = I - D * inv(L)
    Omega = D * Q'
    irfs = L * Q'
    irfs_ortho = L
    return B, Omega, irfs, irfs_ortho
end

# Effect of shock `from` on `outcome` along all paths that go through all
# variables in `var_and` and through none of the variables in `var_not`. Paths
# are enumerated explicitly, making this independent of both methods.
function path_specific_effect(B, Omega, from, outcome, var_and, var_not)
    effect = 0.0
    for mask in 0:(2^(outcome-1)-1)
        nodes = [[i for i in 1:(outcome-1) if (mask >> (i - 1)) & 1 == 1]; outcome]
        all(v -> v in nodes, var_and) || continue
        any(v -> v in nodes, var_not) && continue
        e = Omega[nodes[1], from]
        for t in 1:(length(nodes)-1)
            e *= B[nodes[t+1], nodes[t]]
        end
        effect += e
    end
    return effect
end

@testset "Masking of uninterpretable effects" begin
    B, Omega, irfs, irfs_ortho = simulate_masking_system(4)
    atol = sqrt(eps())

    # sanity check of the simulated system
    @test isapprox((I - B) \ Omega, irfs; atol=atol)
    @test all(isapprox(path_specific_effect(B, Omega, 1, i, Int[], Int[]), irfs[i, 1]; atol=atol) for i in 1:4)

    # Effects of condition `s` using both methods
    function effects(s)
        q = make_condition(s)
        return [
            :BOmega => transmission(1, B, Omega, q; method=:BOmega),
            :irfs => transmission(1, irfs, irfs_ortho, q; method=:irfs),
        ]
    end
    pse(outcome, var_and, var_not) = path_specific_effect(B, Omega, 1, outcome, var_and, var_not)

    @testset "!x1 ($method)" for (method, effect) in effects("!x1")
        @test isapprox(effect[1], 0; atol=atol)
        @test isapprox(effect[2], pse(2, Int[], [1]); atol=atol)
    end

    @testset "x1 ($method)" for (method, effect) in effects("x1")
        @test isnan(effect[1])
        @test isapprox(effect[2], pse(2, [1], Int[]); atol=atol)
    end

    @testset "x2 ($method)" for (method, effect) in effects("x2")
        @test isnan(effect[1])
    end

    @testset "!x2 ($method)" for (method, effect) in effects("!x2")
        @test isapprox(effect[1], pse(1, Int[], [2]); atol=atol)
    end

    @testset "x2 & !x3 ($method)" for (method, effect) in effects("x2 & !x3")
        @test isnan(effect[1])
        @test isnan(effect[2])
        @test isapprox(effect[3], 0; atol=atol)
        @test isapprox(effect[4], pse(4, [2], [3]); atol=atol)
    end

    @testset "!x2 & x3 ($method)" for (method, effect) in effects("!x2 & x3")
        @test isnan(effect[1])
        @test isnan(effect[2])
        @test isnan(effect[3])
        @test isapprox(effect[4], pse(4, [3], [2]); atol=atol)
    end

    # A NOT on the last variable never masks earlier variables (this is the
    # situation of channels such as not_through(model, :x, 0:max_horizon, ...)).
    @testset "!x4 ($method)" for (method, effect) in effects("!x4")
        @test isapprox(effect[1:3], irfs[1:3, 1]; atol=atol)
        @test isapprox(effect[4], 0; atol=atol)
    end

    # Masking is determined by the highest AND variable of any term of the
    # (simplified) condition. x1 | x3 = x1 + x3 - x1 & x3, so x1 to x3 are
    # masked.
    @testset "x1 | x3 ($method)" for (method, effect) in effects("x1 | x3")
        @test all(isnan, effect[1:3])
        @test isapprox(effect[4], pse(4, [1], Int[]) + pse(4, [3], Int[]) - pse(4, [1, 3], Int[]); atol=atol)
    end

    # Both methods must mask the same entries and agree on all others.
    @testset "methods agree: $s" for s in ["!x1", "x1", "x2", "!x2", "x2 & !x3", "!x2 & x3", "!x4", "x1 | x3"]
        (_, e_BOmega), (_, e_irfs) = effects(s)
        @test isnan.(e_BOmega) == isnan.(e_irfs)
        @test isapprox(filter(!isnan, e_BOmega), filter(!isnan, e_irfs); atol=atol)
    end
end
