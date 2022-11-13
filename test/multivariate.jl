using Test
using Distributions
using KernelDensity

import KernelDensity: kernel_dist, default_bandwidth, kde_boundary, kde_range, tabulate

# for D in [Tuple{Normal,Normal}, Tuple{Uniform,Uniform}, Tuple{Logistic,Logistic},
#           (Normal, Normal), (Uniform, Uniform), (Logistic, Logistic)]
#     d = KernelDensity.kernel_dist(D,(0.5,0.5))
#     dx,dy = d
#     @test mean(dx) == 0.0
#     @test mean(dy) == 0.0
#     @test std(dx) ≈ 0.5
#     @test std(dy) ≈ 0.5
# end
 
r = kde_range((-2.0,2.0), 128)
@test step(r) > 0
 
for X in ([0.0], [0.0,0.0], [0.0,0.5], [-0.5:0.1:0.5;])
    w = default_bandwidth(X)
    @test w > 0
    lo, hi = kde_boundary(X, w)
    @test lo < hi
    kr = kde_range((lo, hi), 512)
    @test step(kr) > 0

    for D in (Normal, )
        k1 = tabulate((X, X, X), (kr, kr, kr))
        @test isa(k1, MultivariateKDE{3, R} where R)
        @test size(k1.density) == length.(k1.ranges)
        @test all(>=(0.0), k1.density)
        @test sum(k1.density) * prod(step.(k1.ranges)) ≈ 1.0

        k2 = KernelDensity.conv(k1,kernel_dist(Tuple{D, D, D}, (0.1, 0.1, 0.1)))
        @test isa(k2, MultivariateKDE{3, R} where R)
        @test size(k2.density) == length.(k2.ranges)
        @test all(>=(0.0), k2.density)
        @test sum(k2.density) * prod(step.(k2.ranges)) ≈ 1.0

        k3 = kde((X, X, X); kernel = D)
        @test isa(k3, MultivariateKDE{3, R} where R)
        @test size(k3.density) == length.(k3.ranges)
        @test all(>=(0.0), k3.density)
        @test sum(k3.density) * prod(step.(k3.ranges)) ≈ 1.0

        k4 = kde((X, X, X), (kr, kr, kr); kernel = D)
        @test isa(k4, MultivariateKDE{3, R} where R)
        @test size(k4.density) == length.(k4.ranges)
        @test all(>=(0.0), k4.density)
        @test sum(k4.density) * prod(step.(k4.ranges)) ≈ 1.0

        k5 = kde([X X X]; kernel = D)
        @test isa(k5, MultivariateKDE{3, R} where R)
        @test size(k5.density) == length.(k5.ranges)
        @test all(>=(0.0), k5.density)
        @test sum(k5.density) * prod(step.(k5.ranges)) ≈ 1.0

        k6 = kde([X X X], (kr, kr, kr); kernel = D, weights = fill(1.0 / length(X), length(X)))
        @test k4.density ≈ k6.density
    end
end

for d in 1:5
    data = zeros(10, d)
    k = kde(data, npoints = ntuple(_ -> 16, d))
    @test k isa (MultivariateKDE{d, R} where R)
end

k11 = kde([0.0 0.0 0.0; 1.0 1.0 1.0], (r, r, r), bandwidth = (1, 1, 1), weights = [0, 1])
k12 = kde([1.0 1.0 1.0], (r, r, r), bandwidth = (1, 1, 1))
@test k11.density ≈ k12.density
