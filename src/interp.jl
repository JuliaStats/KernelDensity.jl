import Interpolations: interpolate, scale

mutable struct InterpKDE{K,I} <: AbstractKDE
    kde::K
    itp::I
    InterpKDE{K,I}(kde::K, itp::I) where {K,I} = new{K,I}(kde, itp)
end


function InterpKDE(kde::UnivariateKDE, opts...)
    itp_u = interpolate(kde.density, opts...)
    itp = scale(itp_u, kde.x)
    etp = extrapolate(itp, zero(eltype(kde.density)))
    InterpKDE{typeof(kde),typeof(etp)}(kde, etp)
end
InterpKDE(kde::UnivariateKDE) = InterpKDE(kde, BSpline(Quadratic(Line(OnGrid()))))


function InterpKDE(kde::MultivariateKDE, opts...)
    itp_u = interpolate(kde.density,opts...)
    itp = scale(itp_u, kde.ranges...)
    etp = extrapolate(itp, zero(eltype(kde.density)))
    InterpKDE{typeof(kde),typeof(etp)}(kde, etp)
end
InterpKDE(kde::MultivariateKDE) = InterpKDE(kde, BSpline(Quadratic(Line(OnGrid()))))

pdf(ik::InterpKDE,x::Real...) = ik.itp(x...)
pdf(ik::InterpKDE,xs::AbstractVector) = [ik.itp(x) for x in xs]
function pdf(ik::InterpKDE{K, I}, xs::Vararg{AbstractVector, N}) where
    {N, R, K <: MultivariateKDE{N, R}, I}

    [ik.itp(x...) for x in Iterators.product(xs...)]
end

pdf(k::UnivariateKDE,x) = pdf(InterpKDE(k),x)
pdf(k::MultivariateKDE,x...) = pdf(InterpKDE(k),x...)
