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


function InterpKDE(kde::BivariateKDE, opts...)
    itp_u = interpolate(kde.density,opts...)
    itp = scale(itp_u, kde.x, kde.y)
    etp = extrapolate(itp, zero(eltype(kde.density)))
    InterpKDE{typeof(kde),typeof(etp)}(kde, etp)
end
InterpKDE(kde::BivariateKDE) = InterpKDE(kde::BivariateKDE, BSpline(Quadratic(Line(OnGrid()))))

pdf(ik::InterpKDE,x::Real...) = ik.itp(x...)
pdf(ik::InterpKDE,xs::AbstractVector) = [ik.itp(x) for x in xs]
pdf(ik::InterpKDE,xs::AbstractVector,ys::AbstractVector) = [ik.itp(x,y) for x in xs, y in ys]

pdf(k::UnivariateKDE,x) = pdf(InterpKDE(k),x)
pdf(k::BivariateKDE,x,y) = pdf(InterpKDE(k),x,y)
