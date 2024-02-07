import Interpolations: interpolate, scale

struct InterpKDE{K,I} <: AbstractKDE
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

pdf(k::UnivariateKDE,x) = pdf(InterpKDE(k),x)
pdf(ik::InterpKDE,x::Real...) = ik.itp(x...)
pdf(ik::InterpKDE,xs::AbstractVector) = [ik.itp(x) for x in xs]
Base.broadcasted(::typeof(pdf),k::UnivariateKDE,xs) = InterpKDE(k).itp.(xs)

pdf(k::BivariateKDE,x,y) = pdf(InterpKDE(k),x,y)
pdf(ik::InterpKDE,xs::AbstractVector,ys::AbstractVector) = [ik.itp(x,y) for x in xs, y in ys]
