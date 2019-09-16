module KernelDensity

using DocStringExtensions: TYPEDEF, FIELDS
using StatsBase
using Distributions
using Optim
using Interpolations

import StatsBase: RealVector, RealMatrix
import Distributions: twoπ, pdf
import FFTW: rfft, irfft

export kde, kde_lscv, UnivariateKDE, BivariateKDE, InterpKDE, pdf

abstract type AbstractKDE end

include("univariate.jl")
include("bivariate.jl")
include("interp.jl")

end # module
