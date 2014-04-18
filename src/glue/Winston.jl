using Winston
import Winston: plot, imagesc

function plot(p::FramedPlot, k::UnivariateKDE, args...; kwargs...)
    plot(p, k.x, k.density, args...;  kwargs...)
end

function plot(k::UnivariateKDE, args...; kwargs...)
    plot(Winston.ghf(), k, args...; kwargs...)
end

function imagesc(k::BivariateKDE, clims::Winston.Interval; kwargs...)
    imagesc(extrema(k.x), reverse(extrema(k.y)), k.density', clims; kwargs...)
end
function imagesc(k::BivariateKDE; kwargs...)
    imagesc(k, (0.0,maximum(k.density)); kwargs...)
end
