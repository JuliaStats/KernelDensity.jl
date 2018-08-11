using Test
using Distributions
using KernelDensity

import KernelDensity: kernel_dist, default_bandwidth, kde_boundary, kde_range, tabulate

for D in [Normal,Uniform,Logistic]
    d = kernel_dist(D,0.5)
    @test isa(d,D)
    @test mean(d) == 0.0
    @test std(d) ≈ 0.5
end

r = kde_range((-2.0,2.0), 128)
@test step(r) > 0
r2 = kde_range((0.12698109160784082, 0.9785547869337731), 256)
@test length(r2) == 256

for X in ([0.0], [0.0,0.0], [0.0,0.5], [-0.5:0.1:0.5;])
    w = default_bandwidth(X)
    @test w > 0
    lo, hi = kde_boundary(X,w)
    @test lo < hi
    kr = kde_range((lo,hi), 10)
    @test step(kr) > 0

    for D in (Normal, )
        k1 = tabulate(X,r)
        @test isa(k1,UnivariateKDE)
        @test length(k1.density) == length(k1.x)
        @test all(k1.density .>= 0.0)
        @test sum(k1.density)*step(k1.x) ≈ 1.0

        k2 = KernelDensity.conv(k1,kernel_dist(D,0.1))
        @test isa(k2,UnivariateKDE)
        @test length(k2.density) == length(k2.x)
        @test all(k2.density .>= 0.0)
        @test sum(k2.density)*step(k2.x) ≈ 1.0

        k3 = kde(X;kernel=D)
        @test isa(k3,UnivariateKDE)
        @test length(k3.density) == length(k3.x)
        @test all(k3.density .>= 0.0)
        @test sum(k3.density)*step(k3.x) ≈ 1.0

        k4 = kde(X,r;kernel=D)
        @test isa(k4,UnivariateKDE)
        @test length(k4.density) == length(k4.x)
        @test all(k4.density .>= 0.0)
        @test sum(k4.density)*step(k4.x) ≈ 1.0

        k5 = kde_lscv(X)
        @test isa(k5,UnivariateKDE)
        @test length(k5.density) == length(k5.x)
        @test all(k5.density .>= 0.0)
        @test sum(k5.density)*step(k5.x) ≈ 1.0

        k6 = kde(X,r;kernel=D, weights=fill(1.0/length(X),length(X)))
        @test k4.density ≈ k6.density
    end
end

k11 = kde([0.0, 1.], r, bandwidth=1, weights=[0,1])
k12 = kde([1.], r, bandwidth=1)
@test k11.density ≈ k12.density
