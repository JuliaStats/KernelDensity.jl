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

@testset "pdf method interface" begin
    k = kde([-1., 1.])
    ik = InterpKDE(k)

    @test pdf.(k, [0., 1.]) ≈ [pdf(k, 0.), pdf(k, 1.)] ≈ pdf.(k,[0., 1.]) ≈ [pdf(ik, 0.), pdf(ik, 1.)]
    @test all( pdf.(k, (0., 1.)) .≈ (pdf(k, 0.), pdf(k, 1.)) )
    @test pdf.(k, [0. 1.; 2. -1.]) ≈ [pdf(k, 0.) pdf(k, 1.); pdf(k, 2.) pdf(k, -1.)]

    k2d = kde(([-1., 1.], [0., 1.]))
    ik2d = InterpKDE(k2d)

    @test pdf(k2d, [0.5, 0.1]) ≈ pdf(k2d, 0.5, 0.1) ≈ pdf(ik2d, 0.5, 0.1)
    @test pdf(k2d, [0.5 1.; 0.1 1.]) ≈ [pdf(ik2d, 0.5, 0.1), pdf(ik2d, 1., 1.)]
    @test pdf(k2d, [0.5; 1. ;;; 0.1; 1.]) ≈ [pdf(ik2d, 0.5, 1.) pdf(ik2d, 0.1, 1.)]
end