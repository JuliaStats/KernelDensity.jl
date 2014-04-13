# Store both grid and density for KDE over R2
immutable BivariateKDE
    x::Vector{Float64}
    y::Vector{Float64}
    density::Matrix{Float64}
end


# Algorithm from MASS Chapter 5 for calculating 2D KDE
function kde(x::RealVector, y::RealVector; width::Float64=NaN, resolution::Int=25)
    n = length(x)

    if length(y) != n
        error("x and y must have the same length")
    end

    if isnan(width)
        h1 = bandwidth(x)
        h2 = bandwidth(y)
    else
        h1 = width
        h2 = width
    end

    min_x, max_x = extrema(x)
    min_y, max_y = extrema(y)

    grid_x = [min_x:((max_x - min_x) / (resolution - 1)):max_x]
    grid_y = [min_y:((max_y - min_y) / (resolution - 1)):max_y]

    mx = Array(Float64, resolution, n)
    my = Array(Float64, resolution, n)
    for i in 1:resolution
        for j in 1:n
            mx[i, j] = pdf(Normal(), (grid_x[i] - x[j]) / h1)
            my[i, j] = pdf(Normal(), (grid_y[i] - y[j]) / h2)
        end
    end

    z = A_mul_Bt(mx, my)
    for i in 1:(resolution^2)
        z[i] /= (n * h1 * h2)
    end

    return BivariateKDE(grid_x, grid_y, z)
end
