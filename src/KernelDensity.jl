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

"""n-dimensional convolution"""
function conv(x::AbstractArray{T,N}, w::AbstractArray{T,N}) where {T,N}
	wdim = Int.(ceil.((size(w).-1)./2))
	padding = Iterators.flatten([ (wdim[i],wdim[i]) for i=1:length(wdim) ]) |> collect

	dims = DenseConvDims((size(x)...,1,1),(size(w)...,1,1); padding=padding )
	result = Tracker.conv( reshape(x,(size(x)...,1,1)), reshape(w,(size(w)...,1,1)), dims)
	return dropdims(result, dims = (1+N,2+N))
end

# patches for TrackedReal and Vector{TrackedReal}
conv(x::AbstractArray{TrackedReal{T},N}, w::AbstractArray) where {T,N} = conv(Tracker.collect(x),w)
conv(x::AbstractArray, w::AbstractArray{TrackedReal{T},N}) where {T,N} = conv(x,Tracker.collect(w))
conv(x::AbstractArray{TrackedReal{T},N}, w::AbstractArray{TrackedReal{T},N}) where {T,N} = conv(Tracker.collect(x),Tracker.collect(w))

round(::Type{R}, t::TrackedReal) where {R<:Real} = round(R, t.data)
round(t::TrackedReal, mode::RoundingMode) = round(t.data, mode)
Integer(x::TrackedReal) = Integer(x.data)

include("univariate.jl")
include("bivariate.jl")
include("interp.jl")

end # module
