module KernelDensity

using DocStringExtensions: TYPEDEF, FIELDS
using StatsBase
using Distributions
using Optim
using Interpolations
using Flux
using Flux.Tracker: TrackedReal

import StatsBase: RealVector, RealMatrix
import Distributions: twoπ, pdf
import FFTW: rfft, irfft
import Flux.Tracker: conv
import Base.round

export kde, kde_lscv, UnivariateKDE, BivariateKDE, InterpKDE, pdf

abstract type AbstractKDE end
round(::Type{R}, t::TrackedReal) where {R<:Real} = round(R, t.data)

include("univariate.jl")
include("bivariate.jl")
include("interp.jl")

end # module
