import PyPlot

function PyPlot.plot(k::UnivariateKDE, args...; kwargs...)
    PyPlot.plot(k.x, k.density, args...; kwargs...)
end

function PyPlot.contour(k::BivariateKDE, args...; kwargs...)
    PyPlot.contour(k.x, k.y, k.density', args...; kwargs...)
end
function PyPlot.plot_surface(k::BivariateKDE, args...; kwargs...)
    PyPlot.plot_surface(k.x, k.y, k.density', args...; kwargs...)
end
function PyPlot.surf(k::BivariateKDE, args...; kwargs...)
    PyPlot.surf(k.x, k.y, k.density', args...; kwargs...)
end

function PyPlot.plot_wireframe(k::BivariateKDE, args...; kwargs...)
    PyPlot.plot_wireframe(k.x, k.y, k.density', args...; kwargs...)
end

function PyPlot.imshow(k::BivariateKDE, args...; kwargs...)
    PyPlot.imshow(k.density', args...; origin="lower", extent=[extrema(k.x)...,extrema(k.y)...], kwargs...)
end
