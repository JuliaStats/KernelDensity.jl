import Interpolations: interpolate

type InterpKDE{K,I} <: AbstractKDE
    kde::K
    itp::I
    InterpKDE(kde,itp) = new(kde,itp)
end


function InterpKDE(kde::UnivariateKDE, opts...)
    itp_u = interpolate(kde.density, opts...)
    itp = scale(itp_u, kde.x)
    InterpKDE{typeof(kde),typeof(itp)}(kde, itp)
end
InterpKDE(kde::UnivariateKDE) = InterpKDE(kde, BSpline(Quadratic(Line())), OnGrid())


function InterpKDE(kde::BivariateKDE, opts...)
    itp_u = interpolate(kde.density,opts...)
    itp = scale(itp_u, kde.x, kde.y)
    InterpKDE{typeof(kde),typeof(itp)}(kde, itp)
end
InterpKDE(kde::BivariateKDE) = InterpKDE(kde::BivariateKDE, BSpline(Quadratic(Line())), OnGrid())

pdf(ik::InterpKDE,x::Real...) = ik.itp[x...]
pdf(ik::InterpKDE,xs::AbstractVector) = [ik.itp[x] for x in xs]
pdf(ik::InterpKDE,xs::AbstractVector,ys::AbstractVector) = [ik.itp[x,y] for x in xs, y in ys]

pdf(k::UnivariateKDE,x) = pdf(InterpKDE(k),x)
pdf(k::BivariateKDE,x,y) = pdf(InterpKDE(k),x,y)
