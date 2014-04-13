module KDE

using StatsBase
using Distributions

import Base: conv, FloatRange
import StatsBase: RealVector
import Distributions: twoπ

export kde, UnivariateKDE, BivariateKDE

include("univariate.jl")
include("bivariate.jl")

end # module
