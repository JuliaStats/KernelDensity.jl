

# Silverman's rule of thumb for KDE bandwidth selection
function get_bandwidth(::Silverman, data::AbstractVector{<:Real}, kernelType = Normal, prior = DiscreteUniform(1,length(data)), alpha::Float64 = 0.9)
    # Determine length of data
    ndata = length(data)
    ndata <= 1 && return alpha

    # Calculate width using variance and IQR
    var_width = std(data)
    q25, q75 = quantile(data, (0.25, 0.75))
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
    return alpha * width * ndata^(-0.2), nothing
end

@kwdef struct LSCV <: BandwidthMethod
    nPoints::Int = 2048
    initBandwidth::Float64 = NaN
end

# Select bandwidth using least-squares cross validation, from:
#   Density Estimation for Statistics and Data Analysis
#   B. W. Silverman (1986)
#   sections 3.4.3 (pp. 48-52) and 3.5 (pp. 61-66)
function get_bandwidth(lscv::LSCV, data::AbstractVector{<:Real}, kernelType = Normal, prior = DiscreteUniform(1,length(data)))
    K = lscv.nPoints
    initBandwidth::Float64 = isnan(lscv.initBandwidth) ? get_bandwidth(Silverman(), data)[1] : lscv.initBandwidth
    ndata = length(data)
    lo, hi = extrema(data)
    midpoints = range(lo - 4.0*initBandwidth, hi + 4.0*initBandwidth, K)
    initDen = tabulate(data, midpoints, prior).values

    # the ft here is K/ba*sqrt(2pi) * u(s), it is K times the Yl in Silverman's book
    ft = rfft(initDen)

    ft2 = abs2.(ft)
    c = -twoπ/(step(midpoints)*K)
    hlb, hub = 0.25*initBandwidth, 1.5*initBandwidth

    optimalBandwidth = optimize(hlb, hub) do h
        dist = kernel_dist(kernelType, h)
        ψ = 0.0
        for j in 1:length(ft2)-1
            ks = real(cf(dist, j*c))
            ψ += ft2[j+1]*(ks-2.0)*ks
        end
        ψ*step(midpoints)/K + pdf(dist,0.0)/ndata
    end

    dist = kernel_dist(kernelType, optimalBandwidth)
    for j in 0:length(ft)-1
        ft[j+1] *= cf(dist, j*c)
    end

   convolved = irfft(ft, K)

    # fix rounding error.
    convolved .= max.(0., convolved)

    # Invert the Fourier transform to get the KDE
    optimalBandwidth, DiscretisedPDF(convolved, midpoints)
end