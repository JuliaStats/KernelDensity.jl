const BivariateKDE{Rx <: AbstractRange, Ry <: AbstractRange} =
    MultivariateKDE{2, Tuple{Rx, Ry}}


const BivariateDistribution = NTuple{2, UnivariateDistribution}

function Base.getproperty(k::BivariateKDE, s::Symbol)
    if s === :x
        k.ranges[1]
    elseif s === :y
        k.ranges[2]
    else
        getfield(k, s)
    end
end
