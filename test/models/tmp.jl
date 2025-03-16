using Base: is_file_tracked
k = 3
p = 4
T = 1_000
trend_exponents = [0]

B = 0.2 * randn(k, k*p + length(trend_exponents))
model = simulate(VAR, T, B; trend_exponents=trend_exponents)
fit!(model)
coeffs(model)
B

errors = zeros(k, T)
initial = fill(100, k*p)
model = simulate!(VAR, errors, B; trend_exponents=trend_exponents, initial = initial)
fit!(model)
coeffs(model) - B
maximum(abs, coeffs(model) - B)

model = simulate(VAR, T, B; trend_exponents=trend_exponents)
data = get_intput_data(model)
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

# FIX: This needs fixing
is_stable(model_best)

is_fitted(model)
model = simulate(VAR, T, B; trend_exponents=trend_exponents)
is_fitted(model)
fit!(model)
get_dependent(model)
get_independent(model)
get_intput_data(model)
nobs(model)
residuals(model)
fitted(model)
coeffs(model)


using DataFrames
errors = randn(k, T)
DataFrame(errors', "Y" .* string.(1:k))

