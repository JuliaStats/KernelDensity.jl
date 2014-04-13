# Store both grid and density for KDE over the real line
immutable UnivariateKDE{R<:Range}
    x::R
    density::Vector{Float64}
end

# Silverman's rule of thumb for KDE bandwidth selection
function bandwidth(data::Vector{Float64}, alpha::Float64 = 0.9)
    # Determine length of data
    ndata = length(data)

    # Calculate width using variance and IQR
    var_width = std(data)
    q25, q75 = quantile(data, [0.25, 0.75])
    quantile_width = (q75 - q25) / 1.34

    # Deal with edge cases with 0 IQR or variance
    width = min(var_width, quantile_width)
    if width == 0.0
        if var_width == 0.0
            width = 1.0
        else
            width = var_width
        end
    end

    # Set bandwidth using Silverman's rule of thumb
    return alpha * width * ndata^(-0.2)
end

# construct kernel from bandwidth
kernel_dist(::Type{Normal},w::Real) = Normal(0.0,w)

# Roughly based on:
#   B. W. Silverman (1982) "Algorithm AS 176: Kernel Density Estimation Using
#   the Fast Fourier Transform", Journal of the Royal Statistical
#   Society. Series C (Applied Statistics) , Vol. 31, No. 1, pp. 93-99
#   URL: http://www.jstor.org/stable/2347084
# and:
#   M. C. Jones and H. W. Lotwick (1984) "Remark AS R50: A Remark on Algorithm
#   AS 176. Kernal Density Estimation Using the Fast Fourier Transform",
#   Journal of the Royal Statistical Society. Series C (Applied Statistics) ,
#   Vol. 33, No. 1, pp. 120-122
#   URL: http://www.jstor.org/stable/2347674
# tabulate data for kde

function kde_tab(data::RealVector, npoints::Int=2048)
    # Determine length of data
    ndata = length(data)

    # Double-pad to eliminate cyclic correlation from fft
    lo, hi = extrema(data)
    spread = hi-lo
    lo -= 0.5*spread
    hi += 0.5*spread

    # Set up a grid for discretized data
    grid = zeros(Float64, npoints)

    # Define some more constants
    step = (hi - lo) / npoints
    midpoints = lo:step:hi

    ainc = 1.0 / (ndata*step*step)

    # weighted discetization (cf. Jones and Lotwick)
    for x in data
        k = searchsortedfirst(midpoints,x)
        j = k-1
        if 1 <= j <= npoints
            grid[j] += (midpoints[k]-x)*ainc
            grid[k] += (x-midpoints[j])*ainc
        end
    end

    # returns an un-convolved KDE    
    UnivariateKDE(midpoints, grid)
end


# convolve raw KDE with distribution
function conv(k::UnivariateKDE, d::Distribution)
    # Transform to Fourier basis
    K = length(k.density)
    ft = rfft(k.density)

    # Convolve fft with characteristic function of kernel
    # empirical cf
    #  = \sum_{n=1}^N e^{i*t*X_n} / N
    #  = \sum_{k=0}^K e^{i*t*(a+k*s)} N_k / N
    #  = e^{i*t*a} \sum_{k=0}^K e^{-2pi*i*k*(-t*s*K/2pi)/K} N_k / N
    #  = A * fft(N_k/N)[-t*s*K/2pi + 1]
    c = -twoπ/(step(k.x)*K)
    for j = 1:length(ft)
        ft[j] *= cf(d,(j-1)*c)
    end

    # Invert the Fourier transform to get the KDE
    UnivariateKDE(k.x, irfft(ft, K))
end


function kde(data::RealVector, d::Distribution ; npoints::Int=2048)
    # tabulate data
    k = kde_tab(data, npoints)
    # convolve with distribution
    conv(k,d)
end

function kde(data::RealVector; npoints::Int=2048, width=bandwidth(data), kernel=Normal)    
    width <= 0.0 && error("Bandwidth must be positive")
    d = kernel_dist(kernel,width)
    kde(data,d; npoints=npoints)
end
