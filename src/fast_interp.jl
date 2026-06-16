import FastInterpolations: interp, QuadraticInterp, FillExtrap, AbstractInterpMethod

mutable struct FastInterpKDE{K,I} <: AbstractKDE
    kde::K
    itp::I
    FastInterpKDE{K,I}(kde::K, itp::I) where {K,I} = new{K,I}(kde, itp)
end

function FastInterpKDE(kde::UnivariateKDE; method::AbstractInterpMethod = QuadraticInterp())
    itp = interp(kde.x, kde.density; method, extrap = FillExtrap(0.0))
    FastInterpKDE{typeof(kde), typeof(itp)}(kde, itp)
end

function FastInterpKDE(kde::BivariateKDE; method::AbstractInterpMethod = QuadraticInterp())
    itp = interp((kde.x, kde.y), kde.density; method, extrap = FillExtrap(0.0))
    FastInterpKDE{typeof(kde), typeof(itp)}(kde, itp)
end

pdf(fik::FastInterpKDE, x::Real...)         = fik.itp(x...)
pdf(fik::FastInterpKDE, xs::AbstractVector) = fik.itp(xs)
pdf(fik::FastInterpKDE, xs::AbstractVector, ys::AbstractVector) = [fik.itp(x, y) for x in xs, y in ys]


# One-shot API: a direct query builds no persistent interpolant.
# 1D handles a scalar or a whole vector in a single one-shot call.
function pdf(k::UnivariateKDE, x; method::AbstractInterpMethod = QuadraticInterp())
    return interp(k.x, k.density, x; method = method, extrap = FillExtrap(0.0))
end

# 2D single point: one-shot (a local build, no full-grid interpolant) — the fast path.
function pdf(k::BivariateKDE, x::Real, y::Real; method::AbstractInterpMethod = QuadraticInterp())
    return interp((k.x, k.y), k.density, (x, y); method = method, extrap = FillExtrap(0.0))
end

# 2D grid: build the interpolant once, then evaluate the xs × ys grid (the ND one-shot
# has no efficient outer-product grid form).
function pdf(k::BivariateKDE, xs::AbstractVector, ys::AbstractVector; method::AbstractInterpMethod = QuadraticInterp())
    return pdf(FastInterpKDE(k; method = method), xs, ys)
end