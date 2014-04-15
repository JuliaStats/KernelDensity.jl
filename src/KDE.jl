module KDE

using StatsBase
using Distributions

import Base: conv
import StatsBase: RealVector, RealMatrix
import Distributions: twoπ

export kde, UnivariateKDE, BivariateKDE

include("univariate.jl")
include("bivariate.jl")

end # module
