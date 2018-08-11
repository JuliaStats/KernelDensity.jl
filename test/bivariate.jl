using Test
using Distributions
using KernelDensity

import KernelDensity: kernel_dist, default_bandwidth, kde_boundary, kde_range, tabulate

for D in [Tuple{Normal,Normal}, Tuple{Uniform,Uniform}, Tuple{Logistic,Logistic},
          (Normal, Normal), (Uniform, Uniform), (Logistic, Logistic)]
    d = KernelDensity.kernel_dist(D,(0.5,0.5))
    dx,dy = d
    @test mean(dx) == 0.0
    @test mean(dy) == 0.0
    @test std(dx) ≈ 0.5
    @test std(dy) ≈ 0.5
end

r = kde_range((-2.0,2.0), 128)
@test step(r) > 0

for X in ([0.0], [0.0,0.0], [0.0,0.5], [-0.5:0.1:0.5;])
    w = default_bandwidth(X)
    @test w > 0
    lo, hi = kde_boundary(X,w)
    @test lo < hi
    kr = kde_range((lo,hi), 10)
    @test step(kr) > 0

    for D in (Normal, )
        k1 = tabulate((X,X),(r,r))
        @test isa(k1,BivariateKDE)
        @test size(k1.density) == (length(k1.x), length(k1.y))
        @test all(k1.density .>= 0.0)
        @test sum(k1.density)*step(k1.x)*step(k1.y) ≈ 1.0

        k2 = KernelDensity.conv(k1,kernel_dist(Tuple{D,D}, (0.1,0.1)))
        @test isa(k2,BivariateKDE)
        @test size(k2.density) == (length(k2.x), length(k2.y))
        @test all(k2.density .>= 0.0)
        @test sum(k2.density)*step(k2.x)*step(k2.y) ≈ 1.0

        k3 = kde((X,X);kernel=D)
        @test isa(k3,BivariateKDE)
        @test size(k3.density) == (length(k3.x), length(k3.y))
        @test all(k3.density .>= 0.0)
        @test sum(k3.density)*step(k3.x)*step(k3.y) ≈ 1.0

        k4 = kde((X,X),(r,r);kernel=D)
        @test isa(k4,BivariateKDE)
        @test size(k4.density) == (length(k4.x), length(k4.y))
        @test all(k4.density .>= 0.0)
        @test sum(k4.density)*step(k4.x)*step(k4.y) ≈ 1.0

        k5 = kde([X X];kernel=D)
        @test isa(k5,BivariateKDE)
        @test size(k5.density) == (length(k5.x), length(k5.y))
        @test all(k5.density .>= 0.0)
        @test sum(k5.density)*step(k5.x)*step(k5.y) ≈ 1.0

        k6 = kde([X X],(r,r);kernel=D, weights=fill(1.0/length(X),length(X)))
        @test k4.density ≈ k6.density
    end
end

k11 = kde([0.0 0.0; 1.0 1.0], (r,r), bandwidth=(1,1), weights=[0,1])
k12 = kde([1.0 1.0], (r,r), bandwidth=(1,1))
@test k11.density ≈ k12.density
