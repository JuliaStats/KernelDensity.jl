

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


mutable struct KernelEstimate{VF<:VariateForm, VS<:ValueSupport, KernelType<:Distribution, PriorDist<:UnivariateDistribution{Discrete}, PDF<:DiscretisedPDF} <: AbstractMixtureModel{VF, VS, KernelType}
    const data::Matrix{Float64}
    const prior::PriorDist
    const kernel::KernelType
    const bandwidth::Float64
    precomputedPDF::Union{Nothing, PDF}

    # these guarantee type agreement and nothing more
    function KernelEstimate(data::Matrix{Float64}, prior::PriorDist, kernel::KernelType, bandwidth::Float64, precomputedPDF::PDF) where {
            PriorDist<:UnivariateDistribution{Discrete},
            KernelType <: Distribution,
            PDF <: DiscretisedPDF}
        VF, VS = supertype(KernelType).parameters
        new{VF,VS,KernelType,PriorDist,PDF}(data, prior, kernel, bandwidth, precomputedPDF)
    end
    function KernelEstimate(data::Matrix{Float64}, prior::PriorDist, kernel::KernelType, bandwidth::Float64, ::Nothing) where {
            PriorDist<:UnivariateDistribution{Discrete},
            KernelType <: Distribution}
        VF, VS = supertype(KernelType).parameters
        R = Base.return_types(range,(Float64,Float64,Int))[1] 
        PDF = DiscretisedPDF{size(data)[1],R,eltype(data)}
        new{VF,VS,KernelType,PriorDist,PDF}(data, prior, kernel, bandwidth, nothing)
    end

end

UnivariateKernelEstimate{VS, K, P, PDF} = KernelEstimate{Univariate, VS, K, P, PDF}
MultivariateKernelEstimate{VS, K, P, PDF} = KernelEstimate{Multivariate, VS, K, P, PDF}

abstract type BandwidthMethod end

# default algorithm
struct Silverman<:BandwidthMethod end

# implementing common interface of AbstractMixtureModel
ncomponents(ke::KernelEstimate) = size(ke.data)[2]
component(ke::UnivariateKernelEstimate, k) = ke.kernel - ke.data[1,k]
component(ke::MultivariateKernelEstimate, k) = ke.kernel - ke.data[:,k]
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


# construct kernel from bandwidth
kernel_dist(::Type{Normal}, bandwidth::Real) = Normal(0.0, bandwidth)
kernel_dist(::Type{Uniform}, bandwidth::Real) = Uniform(-√3*bandwidth, √3*bandwidth)

const LocationScale = Union{Laplace, Logistic, SymTriangularDist, Cosine, Epanechnikov}
kernel_dist(::Type{D},w::Real) where {D<:LocationScale} = D(0.0, w/std(D(0.0,1.0)))

# 1D precomputation 

function precompute!(ke::UnivariateKernelEstimate, nPoints::Int = 2048)
    lo, hi = extrema(ke.data)
    midpoints = range(lo - 4.0*ke.bandwidth, hi + 4.0*ke.bandwidth, nPoints)
    ke.precomputedPDF = conv(tabulate(vec(ke.data), midpoints, ke.prior), ke.kernel)
end

function tabulate(data::AbstractVector{<:Real}, midpoints::AbstractRange, prior::UnivariateDistribution{Discrete})
    npoints = length(midpoints)
    s = step(midpoints)

    # Set up a grid for discretized data
    grid = zeros(Float64, npoints)
    ainc = 1.0 / (s*s)

    # weighted discretization (cf. Jones and Lotwick)
    for (i,x) in enumerate(data)
        k = searchsortedfirst(midpoints,x)
        j = k-1
        if 1 <= j <= npoints-1
            grid[j] += (midpoints[k]-x)*ainc*pdf(prior,i)
            grid[k] += (x-midpoints[j])*ainc*pdf(prior,i)
        end
    end

    return DiscretisedPDF(grid, midpoints)
end

function conv(den::DiscretisedPDF{1, R, T}, kernel::UnivariateDistribution) where {T<:Real,R<:AbstractRange}
    # Transform to Fourier basis
    K = length(den.values)
    ft = rfft(den.values)

    # Convolve fft with characteristic function of kernel
    # empirical cf
    #  = \sum_{n=1}^N e^{i*t*X_n} / N
    #  = \sum_{k=0}^K e^{i*t*(a+k*s)} N_k / N
    #  = e^{i*t*a} \sum_{k=0}^K e^{-2pi*i*k*(-t*s*K/2pi)/K} N_k / N
    #  = A * fft(N_k/N)[-t*s*K/2pi + 1]
    c = -twoπ/(step(den.xs)*K)
    for j in 0:length(ft)-1
        ft[j+1] *= cf(kernel,j*c)
    end

    # Invert the Fourier transform to get the KDE
    convolved = irfft(ft, K)

    # fix rounding error.
    convolved .= max.(0., convolved)

    DiscretisedPDF(convolved, den.xs)
end