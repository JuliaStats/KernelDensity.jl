# Store both grid and density for KDE over R2
type BivariateKDE{Rx<:Range,Ry<:Range} <: AbstractKDE
    x::Rx
    y::Ry
    density::Matrix{Float64}
end

function kernel_dist{D<:UnivariateDistribution}(::Type{D},w::(@compat Tuple{Real,Real}))
    kernel_dist(D,w[1]), kernel_dist(D,w[2])
end
function kernel_dist{Dx<:UnivariateDistribution,Dy<:UnivariateDistribution}(::Type{(@compat Tuple{Dx, Dy})},w::(@compat Tuple{Real,Real}))
    kernel_dist(Dx,w[1]), kernel_dist(Dy,w[2])
end

# this function provided for backwards compatibility, though it doesn't have the type restrictions
# to ensure that the given tuple only contains univariate distributions
function kernel_dist(d::(@compat Tuple{DataType, DataType}),w::(@compat Tuple{Real,Real}))
    kernel_dist(d[1],w[1]), kernel_dist(d[2],w[2])
end

# TODO: there are probably better choices.
function default_bandwidth(data::(@compat Tuple{RealVector,RealVector}))
    default_bandwidth(data[1]), default_bandwidth(data[2])
end

# tabulate data for kde
function tabulate(data::(@compat Tuple{RealVector, RealVector}), weights::Weights, midpoints::(@compat Tuple{Range, Range}))
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

function tabulate(data::(@compat Tuple{RealVector, RealVector}), midpoints::(@compat Tuple{Range, Range}))
    tabulate(data, default_weights(data), midpoints)
end

# convolution with product distribution of two univariates distributions
function conv(k::BivariateKDE, dist::(@compat Tuple{UnivariateDistribution,UnivariateDistribution}) )
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

typealias BivariateDistribution @compat(Union{MultivariateDistribution,Tuple{UnivariateDistribution,UnivariateDistribution}})

default_weights(data::(@compat Tuple{RealVector, RealVector})) = UniformWeights(length(data[1]))

function kde(data::(@compat Tuple{RealVector, RealVector}), weights::Weights, midpoints::(@compat Tuple{Range, Range}), dist::BivariateDistribution)
    k = tabulate(data, weights, midpoints)
    conv(k,dist)
end

function kde(data::(@compat Tuple{RealVector, RealVector}), dist::BivariateDistribution;
             boundary::(@compat Tuple{(@compat Tuple{Real,Real}),(@compat Tuple{Real,Real})}) = (kde_boundary(data[1],std(dist[1])),
                                                     kde_boundary(data[2],std(dist[2]))),
             npoints::(@compat Tuple{Int,Int})=(256,256),
             weights::Weights = default_weights(data))

    xmid = kde_range(boundary[1],npoints[1])
    ymid = kde_range(boundary[2],npoints[2])

    kde(data,weights,(xmid,ymid),dist)
end

function kde(data::(@compat Tuple{RealVector, RealVector}), midpoints::(@compat Tuple{Range, Range});
             bandwidth=default_bandwidth(data), kernel=Normal, weights::Weights = default_weights(data))

    dist = kernel_dist(kernel,bandwidth)
    kde(data,weights,midpoints,dist)
end

function kde(data::(@compat Tuple{RealVector, RealVector});
             bandwidth=default_bandwidth(data),
             kernel=Normal,
             boundary::(@compat Tuple{(@compat Tuple{Real,Real}),(@compat Tuple{Real,Real})}) = (kde_boundary(data[1],bandwidth[1]),
                                                     kde_boundary(data[2],bandwidth[2])),
             npoints::(@compat Tuple{Int,Int})=(256,256),
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
