const TrivariateKDE{Rx <: AbstractRange, Ry <: AbstractRange, Rz <: AbstractRange} =
    MultivariateKDE{3, Tuple{Rx, Ry, Rz}}


const TrivariateDistribution = NTuple{3, UnivariateDistribution}

function Base.getproperty(k::TrivariateKDE, s::Symbol)
    if s === :x
        k.ranges[1]
    elseif s === :y
        k.ranges[2]
    elseif s === :z
        k.ranges[3]
    else
        getfield(k, s)
    end
end
