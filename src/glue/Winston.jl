import Winston

function Winston.plot(p::Winston.FramedPlot, k::UnivariateKDE, args...; kwargs...)
    Winston.plot(p, k.x, k.density, args...;  kwargs...)
end

function Winston.plot(k::UnivariateKDE, args...; kwargs...)
    Winston.plot(Winston.ghf(), k, args...; kwargs...)
end

function Winston.imagesc(k::BivariateKDE, clims::Winston.Interval; kwargs...)
    Winston.imagesc(extrema(k.x), reverse(extrema(k.y)), k.density', clims; kwargs...)
end
function Winston.imagesc(k::BivariateKDE; kwargs...)
    Winston.imagesc(k, (0.0,maximum(k.density)); kwargs...)
end
