using Test
using KernelDensity

X = randn(100)
Y = randn(100)

k = kde(X)
@test pdf(k, k.x) ≈ k.density

k = kde((X,Y))
@test pdf(k, k.x, k.y) ≈ k.density

# Try to evaluate the KDE outside the interpolation domain
# The KDE is allowed to be zero, but it should not be greater than the exact solution
k = kde([0.0], bandwidth=1.0)
@test pdf(k, k.x) ≈ k.density
@test pdf(k, -10.0) ≤ pdf(Normal(), -10.0)
@test pdf(k, +10.0) ≤ pdf(Normal(), +10.0)

k = kde(([0.0],[0.0]), bandwidth=(1.0, 1.0))
@test pdf(k, k.x, k.y) ≈ k.density
@test pdf(k, -10.0, 0.0) ≤ pdf(MvNormal(2, 1.0), [-10.0, 0.0])
@test pdf(k, +10.0, 0.0) ≤ pdf(MvNormal(2, 1.0), [+10.0, 0.0])
@test pdf(k, 0.0, -10.0) ≤ pdf(MvNormal(2, 1.0), [0.0, -10.0])
@test pdf(k, 0.0, +10.0) ≤ pdf(MvNormal(2, 1.0), [0.0, +10.0])
