using Base.Test
using KernelDensity
using Grid

srand(123)
X = randn(100)
Y = randn(100)

k = kde(X)
@test_approx_eq pdf(k, k.x) k.density

k = kde((X, Y))
@test_approx_eq pdf(k, k.x, k.y) k.density

# Test BCfill
test_range = [-20:-10; 10:20]
test_len = length(test_range)

k = kde(X)
ik = InterpKDE(k, 0.0, InterpLinear)
@test_approx_eq pdf(ik, k.x) k.density
@test pdf(ik, test_range) == zeros(test_len)

k = kde((X, Y))
ik = InterpKDE(k, 50.0, InterpLinear)
@test_approx_eq pdf(ik, k.x, k.y) k.density
@test_approx_eq pdf(ik, test_range, test_range) fill(50.0, test_len, test_len)
