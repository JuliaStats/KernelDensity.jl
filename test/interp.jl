using Base.Test
using KernelDensity

X = randn(100)
Y = randn(100)

k = kde(X)
@test_approx_eq pdf(k, k.x) k.density

k = kde((X,Y))
@test_approx_eq pdf(k, k.x, k.y) k.density
