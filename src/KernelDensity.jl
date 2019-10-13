module KernelDensity

using DocStringExtensions: TYPEDEF, FIELDS
using StatsBase
using Distributions
using Optim
using Interpolations
using Flux.Tracker: TrackedReal
using Flux

import StatsBase: RealVector, RealMatrix
import Distributions: twoπ, pdf
import FFTW: rfft, irfft
import Flux.Tracker: conv
import Base.round
import Core.Integer

export kde, kde_lscv, UnivariateKDE, BivariateKDE, InterpKDE, pdf

abstract type AbstractKDE end

"""one dimensional convolution using Flux"""
function conv(x::AbstractArray{T,1}, w::AbstractArray{T,1}) where T
	padding = Int(ceil((length(w)-1)/2))
	x,w = reshape(x,(:,1,1)), reshape(w,(:,1,1))

	dims = DenseConvDims(size(x),size(w); padding=(padding,padding))
	conv( x, w, dims)[:,1,1]
end

# patches for TrackedReal and Vector{TrackedReal}
conv(x::AbstractArray{Tracker.TrackedReal{T},1}, w::AbstractArray) where T = conv(Tracker.collect(x),w)
conv(x::AbstractArray, w::AbstractArray{Tracker.TrackedReal{T},1}) where T = conv(x,Tracker.collect(w))
conv(x::AbstractArray{Tracker.TrackedReal{T},1}, w::AbstractArray{Tracker.TrackedReal{T},1}) where T = conv(Tracker.collect(x),Tracker.collect(w))

round(::Type{R}, t::TrackedReal) where {R<:Real} = round(R, t.data)
round(t::TrackedReal, mode::RoundingMode) = round(t.data, mode)
Integer(x::TrackedReal) = Integer(x.data)

include("univariate.jl")
include("bivariate.jl")
include("interp.jl")

end # module
