__precompile__()

module KernelDensity

using StatsBase
using Distributions
using Optim
using Interpolations

import Base: conv
import StatsBase: RealVector, RealMatrix
import Distributions: twoπ, pdf

export kde, kde_lscv, UnivariateKDE, BivariateKDE, InterpKDE, pdf

abstract AbstractKDE

include("univariate.jl")
include("bivariate.jl")
include("interp.jl")

end # module

