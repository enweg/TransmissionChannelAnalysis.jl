using Random

Random.seed!(6150533)

# Recovering the covariance matrix
k = 3
p = 2
T = 10_000_000
trend_exponents = 0:1
B = 0.2 * randn(k, k*p + length(trend_exponents))
Sigma_u = [
    1 0.1 0.1
    0.1 1 0.1
    0.1 0.1 1
]
model = simulate(VAR, T, B, Sigma_u; trend_exponents=trend_exponents)
fit!(model)
@test maximum(abs, cov(model) - Sigma_u) < 1e-2
