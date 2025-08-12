module KernelDensity

using DocStringExtensions: TYPEDEF, FIELDS
using StatsBase
using Distributions
using Interpolations

import Distributions: twoπ, pdf
import FFTW: rfft, irfft

export kde, kde_lscv, UnivariateKDE, BivariateKDE, InterpKDE, pdf

abstract type AbstractKDE end

# define broadcasting behavior
Base.Broadcast.broadcastable(x::AbstractKDE) = Ref(x)

include("univariate.jl")
include("bivariate.jl")
include("interp.jl")

end # module
