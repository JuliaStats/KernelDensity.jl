using Test
using KernelDensity
using Distributions
import FastInterpolations: QuadraticInterp, CubicInterp
import Interpolations

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

@testset "FastInterpKDE" begin
    X = randn(100); Y = randn(100)

    # nodal exactness (1D)
    k = kde(X)
    fik = FastInterpKDE(k)
    @test pdf(fik, k.x) ≈ k.density

    # nodal exactness (2D)
    k2 = kde((X, Y))
    fik2 = FastInterpKDE(k2)
    @test pdf(fik2, k2.x, k2.y) ≈ k2.density

    # tail safety + extrapolation → exactly 0 (1D)
    kc = kde([0.0]; bandwidth = 1.0)
    fkc = FastInterpKDE(kc)
    @test pdf(fkc, kc.x) ≈ kc.density
    @test pdf(fkc, -10.0) == 0.0
    @test pdf(fkc, +10.0) == 0.0
    @test pdf(fkc, -10.0) ≤ pdf(Normal(), -10.0)
    @test pdf(fkc, +10.0) ≤ pdf(Normal(), +10.0)

    # tail safety + extrapolation → exactly 0 (2D)
    kd = kde(([0.0], [0.0]); bandwidth = (1.0, 1.0))
    fkd = FastInterpKDE(kd)
    @test pdf(fkd, kd.x, kd.y) ≈ kd.density
    @test pdf(fkd, -10.0, 0.0) == 0.0
    @test pdf(fkd, 0.0, +10.0) == 0.0

    # method swap stays nodal-exact
    @test pdf(FastInterpKDE(k; method = CubicInterp()), k.x) ≈ k.density
end

@testset "one-shot pdf" begin
    X = randn(100); Y = randn(100)
    k = kde(X); k2 = kde((X, Y))

    # direct queries reproduce the grid nodes
    @test pdf(k, k.x) ≈ k.density                   # 1D vector (one-shot)
    @test pdf(k, 0.37) isa Real                      # 1D scalar (one-shot)
    @test pdf(k2, k2.x, k2.y) ≈ k2.density           # 2D grid (build-once)
    @test pdf(k2, 0.1, 0.2) isa Real                 # 2D scalar (one-shot fast path)

    # extrapolation outside the grid → exactly 0
    @test pdf(k, first(k.x) - 1.0) == 0.0
    @test pdf(k2, first(k2.x) - 1.0, 0.0) == 0.0

    # method keyword threads through
    @test pdf(k, k.x; method = QuadraticInterp()) ≈ k.density
end

@testset "type stability and parity with master (InterpKDE)" begin
    k = kde(randn(80)); k2 = kde((randn(80), randn(80)))

    # 1D: scalar → Float64; vector/range (incl. length-1) → Vector{Float64}
    @test (@inferred pdf(k, 0.5))   isa Float64
    @test (@inferred pdf(k, [0.5])) isa Vector{Float64}
    @test (@inferred pdf(k, k.x))   isa Vector{Float64}

    # 2D: scalar → Float64; vectors/ranges → Matrix{Float64}
    @test (@inferred pdf(k2, 0.1, 0.2))     isa Float64
    @test (@inferred pdf(k2, [0.1], [0.2])) isa Matrix{Float64}
    @test (@inferred pdf(k2, k2.x, k2.y))   isa Matrix{Float64}

    # the same types come out of the InterpKDE path, so the migration is type-transparent
    @test pdf(k, 0.5)             isa Float64
    @test pdf(InterpKDE(k), 0.5)  isa Float64
    @test pdf(k2, k2.x, k2.y)     isa Matrix{Float64}

    # explicit FastInterpKDE wrapper is type-stable too
    @test (@inferred pdf(FastInterpKDE(k), 0.5))       isa Float64
    @test (@inferred pdf(FastInterpKDE(k2), 0.1, 0.2)) isa Float64
end

# The bare pdf(k, …) convenience now routes to the FastInterpolations path, so these
# tests pin the original Interpolations.jl backend explicitly via InterpKDE to guard
# backward compatibility (it is unchanged in src/interp.jl).
@testset "InterpKDE (Interpolations.jl) backward compatibility" begin
    X = randn(100); Y = randn(100)

    # 1D: nodal exactness + tail safety (the original suite's assertions)
    k  = kde(X)
    ik = InterpKDE(k)
    @test pdf(ik, k.x) ≈ k.density

    kc  = kde([0.0]; bandwidth = 1.0)
    ikc = InterpKDE(kc)
    @test pdf(ikc, kc.x) ≈ kc.density
    @test pdf(ikc, -10.0) ≤ pdf(Normal(), -10.0)
    @test pdf(ikc, +10.0) ≤ pdf(Normal(), +10.0)

    # 2D: nodal exactness + tail safety
    k2  = kde((X, Y))
    ik2 = InterpKDE(k2)
    @test pdf(ik2, k2.x, k2.y) ≈ k2.density

    kd  = kde(([0.0], [0.0]); bandwidth = (1.0, 1.0))
    ikd = InterpKDE(kd)
    @test pdf(ikd, kd.x, kd.y) ≈ kd.density
    @test pdf(ikd, -10.0, 0.0) ≤ pdf(MvNormal(2, 1.0), [-10.0, 0.0])
    @test pdf(ikd, 0.0, +10.0) ≤ pdf(MvNormal(2, 1.0), [0.0, +10.0])

    # documented API: extra args are forwarded to Interpolations.interpolate
    bspline = Interpolations.BSpline(Interpolations.Quadratic(Interpolations.Line(Interpolations.OnGrid())))
    @test pdf(InterpKDE(k, bspline), k.x) ≈ k.density

    # return types unchanged: scalar / vector / matrix
    @test pdf(ik, 0.5)         isa Float64
    @test pdf(ik, k.x)         isa Vector{Float64}
    @test pdf(ik2, 0.1, 0.2)   isa Float64
    @test pdf(ik2, k2.x, k2.y) isa Matrix{Float64}
end
