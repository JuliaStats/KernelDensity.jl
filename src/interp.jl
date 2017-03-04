import Interpolations: interpolate, ExtrapDimSpec

type InterpKDE{K,I} <: AbstractKDE
    kde::K
    itp::I
    InterpKDE(kde,itp) = new(kde,itp)
end


function InterpKDE(kde::UnivariateKDE, extrap::Union{ExtrapDimSpec, Number}, opts...)
    itp_u = interpolate(kde.density, opts...)
    itp_u = extrapolate(itp_u, extrap)
    itp = scale(itp_u, kde.x)
    InterpKDE{typeof(kde),typeof(itp)}(kde, itp)
end
InterpKDE(kde::UnivariateKDE) = InterpKDE(kde, NaN, BSpline(Quadratic(Line())), OnGrid())


function InterpKDE(kde::BivariateKDE, extrap::Union{ExtrapDimSpec, Number}, opts...)
    itp_u = interpolate(kde.density,opts...)
    itp_u = extrapolate(itp_u, extrap)
    itp = scale(itp_u, kde.x, kde.y)
    InterpKDE{typeof(kde),typeof(itp)}(kde, itp)
end
InterpKDE(kde::BivariateKDE) = InterpKDE(kde::BivariateKDE, NaN, BSpline(Quadratic(Line())), OnGrid())

pdf(ik::InterpKDE,x::Real...) = ik.itp[x...]
pdf(ik::InterpKDE,xs::AbstractVector) = [ik.itp[x] for x in xs]
pdf(ik::InterpKDE,xs::AbstractVector,ys::AbstractVector) = [ik.itp[x,y] for x in xs, y in ys]

pdf(k::UnivariateKDE,x) = pdf(InterpKDE(k),x)
pdf(k::BivariateKDE,x,y) = pdf(InterpKDE(k),x,y)
