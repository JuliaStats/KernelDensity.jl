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

# interface implementation
# it should be consistent with Distributions.pdf

pdf(k::UnivariateKDE,x) = pdf(InterpKDE(k),x)
Base.broadcasted(::typeof(pdf),k::UnivariateKDE,xs) = Base.broadcasted(InterpKDE(k).itp, xs)
pdf(ik::InterpKDE,xs::AbstractVector) = pdf.(ik, xs)

pdf(k::BivariateKDE,x,y) = pdf(InterpKDE(k),x,y)
pdf(ik::InterpKDE,xs::AbstractVector,ys::AbstractVector) = ik.itp.(xs,ys')
pdf(k::BivariateKDE, M) = pdf(InterpKDE(k),M)
pdf(ik::InterpKDE, M::AbstractArray{<:Real, 1}) = ik.itp(M[1],M[2])
pdf(ik::InterpKDE, M::AbstractArray{<:Real, N}) where N = pdf.(ik,eachslice(M, dims=ntuple(i->i+1, N-1)) )
