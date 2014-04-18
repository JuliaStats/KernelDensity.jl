module KernelDensity

using StatsBase
using Distributions

import Base: conv
import StatsBase: RealVector, RealMatrix
import Distributions: twoπ

export kde, UnivariateKDE, BivariateKDE

include("univariate.jl")
include("bivariate.jl")

macro glue(pkg)
    path = joinpath(dirname(Base.source_path(nothing)),"glue",string(pkg,".jl"))
    init = symbol(string(pkg,"_init"))
    quote
        $(esc(init))() = include($path)
        isdefined(Main,$(QuoteNode(pkg))) && $(esc(init))()
    end
end

@glue Winston

end # module
