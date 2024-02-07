"""
$(TYPEDEF)

Store both grid and density for KDE over ``ℝ²``.

Reading the fields directly is part of the API, and

```julia
sum(density) * step(x) ≈ 1
```

# Fields

$(FIELDS)
"""
mutable struct UnivariateKDE{R<:AbstractRange} <: AbstractKDE
    "Gridpoints for evaluating the density."
    x::R
    "Kernel density at corresponding gridpoints `x`."
    density::Vector{Float64}
end

# construct kernel from bandwidth
kernel_dist(::Type{Normal},w::Real) = Normal(0.0,w)
kernel_dist(::Type{Uniform},w::Real) = (s = 1.7320508075688772*w; Uniform(-s,s))

const LocationScale = Union{Laplace,Logistic,SymTriangularDist}
kernel_dist(::Type{D},w::Real) where {D} = (s = w/std(D(0.0,1.0)); D(0.0,s))


# Silverman's rule of thumb for KDE bandwidth selection
function default_bandwidth(data::AbstractVector{<:Real}, alpha::Float64 = 0.9)
    # Determine length of data
    ndata = length(data)
    ndata <= 1 && return alpha

    # Calculate width using variance and IQR
    var_width = std(data)
    q25, q75 = quantile(data, [0.25, 0.75])
    quantile_width = (q75 - q25) / 1.34

    # Deal with edge cases with 0 IQR or variance
    width = min(var_width, quantile_width)
    if width == 0.0
        if var_width == 0.0
            width = 1.0
        else
            width = var_width
        end
    end

    # Set bandwidth using Silverman's rule of thumb
    return alpha * width * ndata^(-0.2)
end

function default_weights(data::AbstractVector{<:Real})
    UniformWeights(length(data))
end


# Roughly based on:
#   B. W. Silverman (1982) "Algorithm AS 176: Kernel Density Estimation Using
#   the Fast Fourier Transform", Journal of the Royal Statistical
#   Society. Series C (Applied Statistics) , Vol. 31, No. 1, pp. 93-99
#   URL: http://www.jstor.org/stable/2347084
# and:
#   M. C. Jones and H. W. Lotwick (1984) "Remark AS R50: A Remark on Algorithm
#   AS 176. Kernal Density Estimation Using the Fast Fourier Transform",
#   Journal of the Royal Statistical Society. Series C (Applied Statistics) ,
#   Vol. 33, No. 1, pp. 120-122
#   URL: http://www.jstor.org/stable/2347674

# default kde range
# Should extend enough beyond the data range to avoid cyclic correlation from the FFT
function kde_boundary(data::AbstractVector{<:Real}, bandwidth::Real)
    lo, hi = extrema(data)
    lo - 4.0*bandwidth, hi + 4.0*bandwidth
end

# convert boundary and npoints to Range object
function kde_range(boundary::Tuple{Real,Real}, npoints::Int)
    lo, hi = boundary
    lo < hi || error("boundary (a,b) must have a < b")

    range(lo, stop=hi, length=npoints)
end

struct UniformWeights{N} end

UniformWeights(n) = UniformWeights{n}()

Base.sum(x::UniformWeights) = 1.
Base.getindex(x::UniformWeights{N}, i) where {N} = 1/N

const Weights = Union{UniformWeights, AbstractVector{<:Real}, StatsBase.Weights}


# tabulate data for kde
function tabulate(data::AbstractVector{<:Real}, midpoints::R, weights::Weights=default_weights(data)) where R<:AbstractRange
    npoints = length(midpoints)
    s = step(midpoints)

    # Set up a grid for discretized data
    grid = zeros(Float64, npoints)
    ainc = 1.0 / (sum(weights)*s*s)

    # weighted discretization (cf. Jones and Lotwick)
    for (i,x) in enumerate(data)
        k = searchsortedfirst(midpoints,x)
        j = k-1
        if 1 <= j <= npoints-1
            grid[j] += (midpoints[k]-x)*ainc*weights[i]
            grid[k] += (x-midpoints[j])*ainc*weights[i]
        end
    end

    # returns an un-convolved KDE
    UnivariateKDE(midpoints, grid)
end

# convolve raw KDE with kernel
# TODO: use in-place fft
function conv(k::UnivariateKDE, dist::UnivariateDistribution)
    # Transform to Fourier basis
    K = length(k.density)
    ft = rfft(k.density)

    # Convolve fft with characteristic function of kernel
    # empirical cf
    #  = \sum_{n=1}^N e^{i*t*X_n} / N
    #  = \sum_{k=0}^K e^{i*t*(a+k*s)} N_k / N
    #  = e^{i*t*a} \sum_{k=0}^K e^{-2pi*i*k*(-t*s*K/2pi)/K} N_k / N
    #  = A * fft(N_k/N)[-t*s*K/2pi + 1]
    c = -twoπ/(step(k.x)*K)
    for j = 0:length(ft)-1
        ft[j+1] *= cf(dist,j*c)
    end

    dens = irfft(ft, K)
    # fix rounding error.
    for i = 1:K
        dens[i] = max(0.0,dens[i])
    end

    # Invert the Fourier transform to get the KDE
    UnivariateKDE(k.x, dens)
end

# main kde interface methods

"""

	kde(data; kwargs...)
	kde((xdata, ydata); kwargs...)
	
Kernel density estimation method. Returns 1D or 2D KDE object. The grid used and the values of the estimated density can be obtained from fields `.x` and `.density` respectively. To obtain kde values at points different than the initial grid use the `pdf` method.

The keyword arguments are
* `boundary`: the lower and upper limits of the kde, tuple in 1D case, tuple of tuples in 2D case, 
* `npoints`: the number of interpolation points to use,
* `kernel = Normal`: the distributional family from [Distributions.jl](https://github.com/JuliaStats/Distributions.jl),
* `bandwidth`: the bandwidth of the kernel; default is calculated using Silverman's rule. 
"""
function kde(data::AbstractVector{<:Real}, weights::Weights, midpoints::R, dist::UnivariateDistribution) where R<:AbstractRange
    k = tabulate(data, midpoints, weights)
    conv(k,dist)
end

function kde(data::AbstractVector{<:Real}, dist::UnivariateDistribution;
             boundary::Tuple{Real,Real}=kde_boundary(data,std(dist)), npoints::Int=2048, weights=default_weights(data))

    midpoints = kde_range(boundary,npoints)
    kde(data,weights,midpoints,dist)
end

function kde(data::AbstractVector{<:Real}, midpoints::R;
             bandwidth=default_bandwidth(data), kernel=Normal, weights=default_weights(data)) where R<:AbstractRange
    bandwidth > 0.0 || error("Bandwidth must be positive")
    dist = kernel_dist(kernel,bandwidth)
    kde(data,weights,midpoints,dist)
end

function kde(data::AbstractVector{<:Real}; bandwidth=default_bandwidth(data), kernel=Normal,
             npoints::Int=2048, boundary::Tuple{Real,Real}=kde_boundary(data,bandwidth), weights=default_weights(data))
    bandwidth > 0.0 || error("Bandwidth must be positive")
    dist = kernel_dist(kernel,bandwidth)
    kde(data,dist;boundary=boundary,npoints=npoints,weights=weights)
end

"""
    optimize(f, x_lower, x_upper; iterations=1000, rel_tol=nothing, abs_tol=nothing)

Minimize the function `f` in the interval `x_lower..x_upper`, using the
[golden-section search](https://en.wikipedia.org/wiki/Golden-section_search).
Return an approximate minimum `x̃` or error if such approximate minimum cannot be found.

This algorithm assumes that `-f` is unimodal on the interval `x_lower..x_upper`,
that is to say, there exists a unique `x` in `x_lower..x_upper` such that `f` is
decreasing on `x_lower..x` and increasing on `x..x_upper`.

`rel_tol` and `abs_tol` determine the relative and absolute tolerance, that is
to say, the returned value `x̃` should differ from the actual minimum `x` at most
`abs_tol + rel_tol * abs(x̃)`.
If not manually specified, `rel_tol` and `abs_tol` default to `sqrt(eps(T))` and
`eps(T)` respectively, where `T` is the floating point type of `x_lower` and `x_upper`.

`iterations` determines the maximum number of iterations allowed before convergence.

This is a private, unexported function, used internally to select the optimal bandwidth
automatically.
"""
function optimize(f, x_lower, x_upper; iterations=1000, rel_tol=nothing, abs_tol=nothing)

    if x_lower > x_upper
        error("x_lower must be less than x_upper")
    end

    T = promote_type(typeof(x_lower/1), typeof(x_upper/1))
    rtol = something(rel_tol, sqrt(eps(T)))
    atol = something(abs_tol, eps(T))
    
    function midpoint_and_convergence(lower, upper)
        midpoint = (lower + upper) / 2
        tol = atol + rtol * midpoint
        midpoint, (upper - lower) <= 2tol
    end
    
    invphi::T = 0.5 * (sqrt(5) - 1)
    invphisq::T = 0.5 * (3 - sqrt(5))
    
    a::T, b::T = x_lower, x_upper
    h = b - a
    c = a + invphisq * h
    d = a + invphi * h
    
    fc, fd = f(c), f(d)
    
    for _ in 1:1000
        h *= invphi
        if fc < fd
            m, converged = midpoint_and_convergence(a, d)
            converged && return m
            b = d
            d, fd = c, fc
            c = a + invphisq * h
            fc = f(c)
        else
            m, converged = midpoint_and_convergence(c, b)
            converged && return m
            a = c
            c, fc = d, fd
            d = a + invphi * h
            fd = f(d)
        end
    end

    error("Reached maximum number of iterations without convergence.")
end

# Select bandwidth using least-squares cross validation, from:
#   Density Estimation for Statistics and Data Analysis
#   B. W. Silverman (1986)
#   sections 3.4.3 (pp. 48-52) and 3.5 (pp. 61-66)

function kde_lscv(data::AbstractVector{<:Real}, midpoints::R;
                  kernel=Normal,
                  bandwidth_range::Tuple{Real,Real}=(h=default_bandwidth(data); (0.25*h,1.5*h)),
                  weights=default_weights(data)) where R<:AbstractRange

    ndata = length(data)
    k = tabulate(data, midpoints, weights)

    # the ft here is K/ba*sqrt(2pi) * u(s), it is K times the Yl in Silverman's book
    K = length(k.density)
    ft = rfft(k.density)

    ft2 = abs2.(ft)
    c = -twoπ/(step(k.x)*K)
    hlb, hub = bandwidth_range

    minimizer = optimize(hlb, hub) do h
        dist = kernel_dist(kernel, h)
        ψ = 0.0
        for j = 1:length(ft2)-1
            ks = real(cf(dist, j*c))
            ψ += ft2[j+1]*(ks-2.0)*ks
        end
        ψ*step(k.x)/K + pdf(dist,0.0)/ndata
    end

    dist = kernel_dist(kernel, minimizer)
    for j = 0:length(ft)-1
        ft[j+1] *= cf(dist, j*c)
    end

    dens = irfft(ft, K)
    # fix rounding error.
    for i = 1:K
        dens[i] = max(0.0,dens[i])
    end

    # Invert the Fourier transform to get the KDE
    UnivariateKDE(k.x, dens)
end

function kde_lscv(data::AbstractVector{<:Real};
                  boundary::Tuple{Real,Real}=kde_boundary(data,default_bandwidth(data)),
                  npoints::Int=2048,
                  kernel=Normal,
                  bandwidth_range::Tuple{Real,Real}=(h=default_bandwidth(data); (0.25*h,1.5*h)),
                  weights::Weights = default_weights(data))

    midpoints = kde_range(boundary,npoints)
    kde_lscv(data,midpoints; kernel=kernel, bandwidth_range=bandwidth_range, weights=weights)
end
