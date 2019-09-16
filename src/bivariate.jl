"""
$(TYPEDEF)

Store both grid and density for KDE over the real line.

Reading the fields directly is part of the API, and

```julia
sum(density) * step(x) * step(y) ≈ 1
```

# Fields

$(FIELDS)
"""
mutable struct BivariateKDE{Rx<:AbstractRange,Ry<:AbstractRange} <: AbstractKDE
    "First coordinate of gridpoints for evaluating the density."
    x::Rx
    "Second coordinate of gridpoints for evaluating the density."
    y::Ry
    "Kernel density at corresponding gridpoints `Tuple.(x, permutedims(y))`."
    density::Matrix{Float64}
end

function kernel_dist(::Type{D},w::Tuple{Real,Real}) where D<:UnivariateDistribution
    kernel_dist(D,w[1]), kernel_dist(D,w[2])
end
function kernel_dist(::Type{Tuple{Dx, Dy}},w::Tuple{Real,Real}) where {Dx<:UnivariateDistribution,Dy<:UnivariateDistribution}
    kernel_dist(Dx,w[1]), kernel_dist(Dy,w[2])
end

const DataTypeOrUnionAll = Union{DataType, UnionAll}

# this function provided for backwards compatibility, though it doesn't have the type restrictions
# to ensure that the given tuple only contains univariate distributions
function kernel_dist(d::Tuple{DataTypeOrUnionAll, DataTypeOrUnionAll}, w::Tuple{Real,Real})
    kernel_dist(d[1],w[1]), kernel_dist(d[2],w[2])
end

# TODO: there are probably better choices.
function default_bandwidth(data::Tuple{RealVector,RealVector})
    default_bandwidth(data[1]), default_bandwidth(data[2])
end

# tabulate data for kde
function tabulate(data::Tuple{RealVector, RealVector}, midpoints::Tuple{Rx, Ry},
        weights::Weights = default_weights(data)) where {Rx<:AbstractRange,Ry<:AbstractRange}
    xdata, ydata = data
    ndata = length(xdata)
    length(ydata) == ndata || error("data vectors must be of same length")

    xmid, ymid = midpoints
    nx, ny = length(xmid), length(ymid)
    sx, sy = step(xmid), step(ymid)

    # Set up a grid for discretized data
    grid = zeros(Float64, nx, ny)
    ainc = 1.0 / (sum(weights)*(sx*sy)^2)

    # weighted discretization (cf. Jones and Lotwick)
    for i in 1:length(xdata)
        x = xdata[i]
        y = ydata[i]
        kx, ky = searchsortedfirst(xmid,x), searchsortedfirst(ymid,y)
        jx, jy = kx-1, ky-1
        if 1 <= jx <= nx-1 && 1 <= jy <= ny-1
            grid[jx,jy] += (xmid[kx]-x)*(ymid[ky]-y)*ainc*weights[i]
            grid[kx,jy] += (x-xmid[jx])*(ymid[ky]-y)*ainc*weights[i]
            grid[jx,ky] += (xmid[kx]-x)*(y-ymid[jy])*ainc*weights[i]
            grid[kx,ky] += (x-xmid[jx])*(y-ymid[jy])*ainc*weights[i]
        end
    end

    # returns an un-convolved KDE
    BivariateKDE(xmid, ymid, grid)
end

# convolution with product distribution of two univariates distributions
function conv(k::BivariateKDE, dist::Tuple{UnivariateDistribution,UnivariateDistribution})
    # Transform to Fourier basis
    Kx, Ky = size(k.density)
    ft = rfft(k.density)

    distx, disty = dist

    # Convolve fft with characteristic function of kernel
    cx = -twoπ/(step(k.x)*Kx)
    cy = -twoπ/(step(k.y)*Ky)
    for j = 0:size(ft,2)-1
        for i = 0:size(ft,1)-1
            ft[i+1,j+1] *= cf(distx,i*cx)*cf(disty,min(j,Ky-j)*cy)
        end
    end
    dens = irfft(ft, Kx)

    for i = 1:length(dens)
        dens[i] = max(0.0,dens[i])
    end

    # Invert the Fourier transform to get the KDE
    BivariateKDE(k.x, k.y, dens)
end

const BivariateDistribution = Union{MultivariateDistribution,Tuple{UnivariateDistribution,UnivariateDistribution}}

default_weights(data::Tuple{RealVector, RealVector}) = UniformWeights(length(data[1]))

function kde(data::Tuple{RealVector, RealVector}, weights::Weights, midpoints::Tuple{Rx, Ry},
        dist::BivariateDistribution) where {Rx<:AbstractRange,Ry<:AbstractRange}
    k = tabulate(data, midpoints, weights)
    conv(k,dist)
end

function kde(data::Tuple{RealVector, RealVector}, dist::BivariateDistribution;
             boundary::Tuple{Tuple{Real,Real}, Tuple{Real,Real}} = (kde_boundary(data[1],std(dist[1])),
                                                     kde_boundary(data[2],std(dist[2]))),
             npoints::Tuple{Int,Int}=(256,256),
             weights::Weights = default_weights(data))

    xmid = kde_range(boundary[1],npoints[1])
    ymid = kde_range(boundary[2],npoints[2])

    kde(data,weights,(xmid,ymid),dist)
end

function kde(data::Tuple{RealVector, RealVector}, midpoints::Tuple{Rx, Ry};
             bandwidth=default_bandwidth(data), kernel=Normal,
             weights::Weights = default_weights(data)) where {Rx<:AbstractRange,Ry<:AbstractRange}

    dist = kernel_dist(kernel,bandwidth)
    kde(data,weights,midpoints,dist)
end

function kde(data::Tuple{RealVector, RealVector};
             bandwidth=default_bandwidth(data),
             kernel=Normal,
             boundary::Tuple{Tuple{Real,Real}, Tuple{Real,Real}} = (kde_boundary(data[1],bandwidth[1]),
                                                     kde_boundary(data[2],bandwidth[2])),
             npoints::Tuple{Int,Int}=(256,256),
             weights::Weights = default_weights(data))

    dist = kernel_dist(kernel,bandwidth)
    xmid = kde_range(boundary[1],npoints[1])
    ymid = kde_range(boundary[2],npoints[2])

    kde(data,weights,(xmid,ymid),dist)
end

# matrix data
function kde(data::RealMatrix,args...;kwargs...)
    size(data,2) == 2 || error("Can only construct KDE from matrices with 2 columns.")
    kde((data[:,1],data[:,2]),args...;kwargs...)
end
