module KernelDensity

using DocStringExtensions: TYPEDEF, FIELDS
using StatsBase
using Distributions
using Interpolations

import IrrationalConstants: twoπ
import Base: getproperty, propertynames

import Distributions: pdf
import Distributions: ncomponents, component, probs

import FFTW: rfft, irfft

export DiscretisedPDF, KernelEstimate, BandwidthMethod, Silverman, LSCV
export kernel_estimate, precompute!

export kde, kde_lscv, UnivariateKDE, BivariateKDE, InterpKDE, pdf

abstract type AbstractKDE end

Base.broadcastable(x::AbstractKDE) = Ref(x)

include("univariate.jl")
include("bivariate.jl")
include("interp.jl")

include("initialisation.jl")
include("bandwidth_selection.jl")
include("feature_computation.jl")

end
