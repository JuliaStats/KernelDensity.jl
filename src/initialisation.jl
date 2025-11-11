

struct DiscretisedPDF{N,R,T}
    values::Array{T,N}
    labels::NTuple{N,R}
    DiscretisedPDF(values::Array{T,N}, labels::NTuple{N,R}) where {N,T<:Real,R<:AbstractRange} = new{N,R,T}(values, labels)
end

DiscretisedPDF(values::Vector, label::AbstractRange) = DiscretisedPDF(values, (label,))

# provides d.xs, d.ys, d.zs for convenience
function Base.getproperty(d::DiscretisedPDF, prop::Symbol)
    if prop == :xs
        return d.labels[1]
    elseif prop == :ys
        return d.labels[2]
    elseif prop == :zs
        return d.labels[3]
    else getfield(d, prop) 
    end
end
Base.propertynames(::DiscretisedPDF) = (:values, :labels, :xs, :ys, :zs)

if VERSION >= v"1.8" # for backwards compatibility no const fields in < v1.8
    include("ke-1.8.jl")
else
    include("ke-1.6.jl")
end


UnivariateKernelEstimate{VS, K, P, PDF} = KernelEstimate{Univariate, VS, K, P, PDF}
MultivariateKernelEstimate{VS, K, P, PDF} = KernelEstimate{Multivariate, VS, K, P, PDF}

# KernelEstimate is a scalar
Base.broadcastable(ke::KernelEstimate) = Ref(ke)

# It is possible to add linear transformations a*ke + b

abstract type BandwidthMethod end

# default algorithm
struct Silverman<:BandwidthMethod end

# implementing common interface of AbstractMixtureModel
ncomponents(ke::KernelEstimate) = size(ke.data)[2]
component(ke::UnivariateKernelEstimate, k) = ke.kernel + ke.data[1,k]
component(ke::MultivariateKernelEstimate, k) = ke.kernel + ke.data[:,k]
probs(ke::KernelEstimate) = probs(ke.prior)


# creating KernelEstimate instance

# make kernel density given bandwidth
function kernel_estimate(data::Vector{<:Real}, bandwidth::Real, kernelType = Normal, prior::UnivariateDistribution{Discrete} = DiscreteUniform(1,length(data)))
    kernel = kernel_dist(kernelType, bandwidth)
    KernelEstimate(reshape(data, 1, length(data)), prior, kernel, Float64(bandwidth), nothing)
end

# find bandwidth, then make kernel density
function kernel_estimate(data::Vector{<:Real}, method::BandwidthMethod = Silverman(), kernelType = Normal, prior::UnivariateDistribution{Discrete} = DiscreteUniform(1,length(data)))
    bandwidth, kde = get_bandwidth(method, data, kernelType, prior)
    kernel = kernel_dist(kernelType, bandwidth)
    KernelEstimate(reshape(data,1,length(data)), prior, kernel, bandwidth, kde)
end

# Can add kernel_estimate which takes prior as a vector. 

# construct kernel from bandwidth
kernel_dist(::Type{Normal}, bandwidth::Real) = Normal(0.0, bandwidth)
kernel_dist(::Type{Uniform}, bandwidth::Real) = Uniform(-√3*bandwidth, √3*bandwidth)

const LocationScale = Union{Laplace, Logistic, SymTriangularDist, Cosine, Epanechnikov}
kernel_dist(::Type{D},w::Real) where {D<:LocationScale} = D(0.0, w/std(D(0.0,1.0)))
