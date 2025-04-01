k = 3
p = 2
T = 1_000
trend_exponents = [0]

B = 0.2 * randn(k, k*p + length(trend_exponents))
model = simulate(VAR, T, B; trend_exponents=trend_exponents)
coeffs(model)
cov(model)
fitted(model)
residuals(model)
aic(model)
fit!(model)
coeffs(model)
B

errors = zeros(k, T)
initial = fill(100, k*p)
model = simulate!(VAR, errors, B; trend_exponents=trend_exponents, initial = initial)
fit!(model)
coeffs(model) - B
maximum(abs, coeffs(model) - B)
cov(model)

model = simulate(VAR, T, B; trend_exponents=trend_exponents)
data = get_input_data(model)
model = VAR(data, 10; trend_exponents=trend_exponents)
model_best, ic_table = fit_and_select!(model, aic)
model_best.p
ic_table
model_best, ic_table = fit_and_select!(model, bic)
model_best.p
ic_table
model_best, ic_table = fit_and_select!(model, hqc)
model_best.p
ic_table
model_best, ic_table = fit_and_select!(model, sic)
model_best.p
ic_table

is_stable(model_best)
make_companion_matrix(model_best)

model = simulate(VAR, T, B; trend_exponents=trend_exponents)
is_fitted(model)
fit!(model)
is_fitted(model)

get_dependent(model)
get_independent(model)
get_input_data(model)
nobs(model)
residuals(model)
fitted(model)
coeffs(model)

k = 3
p = 2
T = 1_000
trend_exponents = 0:1
B = 0.2 * randn(k, k*p + length(trend_exponents))
model = simulate(VAR, T, B; trend_exponents=trend_exponents)

k = 3
p = 2
T = 10_000_000
trend_exponents = 0:1
B = 0.2 * randn(k, k*p + length(trend_exponents))
Sigma_u = randn(k, k)
Sigma_u = (Sigma_u * Sigma_u') / 2
model = simulate(VAR, T, B, Sigma_u; trend_exponents=trend_exponents)
fit!(model)
cov(model) - Sigma_u
