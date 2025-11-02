
# 1D pdf precomputation 

function precompute!(ke::UnivariateKernelEstimate, nPoints::Integer = 2048, 
        boundary::Tuple{<:Real,<:Real} =((lo, hi) = extrema(ke.data); (lo - 4.0*ke.bandwidth, hi + 4.0*ke.bandwidth)))
    
    # find the element type of range in ke.precomputedPDF
    T = eltype(typeof(ke).parameters[5].parameters[2])
    midpoints = range(T(boundary[1]), T(boundary[2]), Int(nPoints))
    ke.precomputedPDF = conv(tabulate(vec(ke.data), midpoints, ke.prior), ke.kernel)
end

function tabulate(data::AbstractVector{<:Real}, midpoints::AbstractRange, prior::UnivariateDistribution{Discrete})
    npoints = length(midpoints)
    s = step(midpoints)

    # Set up a grid for discretized data
    grid = zeros(Float64, npoints)
    ainc = 1.0 / (s*s)

    # weighted discretization (cf. Jones and Lotwick)
    for (i,x) in enumerate(data)
        k = searchsortedfirst(midpoints,x)
        j = k-1
        if 1 <= j <= npoints-1
            grid[j] += (midpoints[k]-x)*ainc*pdf(prior,i)
            grid[k] += (x-midpoints[j])*ainc*pdf(prior,i)
        end
    end

    return DiscretisedPDF(grid, midpoints)
end

function conv(den::DiscretisedPDF{1, R, T}, kernel::UnivariateDistribution) where {T<:Real,R<:AbstractRange}
    # Transform to Fourier basis
    K = length(den.values)
    ft = rfft(den.values)

    # Convolve fft with characteristic function of kernel
    # empirical cf
    #  = \sum_{n=1}^N e^{i*t*X_n} / N
    #  = \sum_{k=0}^K e^{i*t*(a+k*s)} N_k / N
    #  = e^{i*t*a} \sum_{k=0}^K e^{-2pi*i*k*(-t*s*K/2pi)/K} N_k / N
    #  = A * fft(N_k/N)[-t*s*K/2pi + 1]
    c = -twoπ/(step(den.xs)*K)
    for j in 0:length(ft)-1
        ft[j+1] *= cf(kernel,j*c)
    end

    # Invert the Fourier transform to get the KDE
    convolved = irfft(ft, K)

    # fix rounding error.
    convolved .= max.(0., convolved)

    DiscretisedPDF(convolved, den.xs)
end

# pdf methods

# pdf(ke, x) is implemented in Distributions.jl 

function pdf(ke::UnivariateKernelEstimate, x::Real, method::Symbol)
    if method == :precomputed
        den = ke.precomputedPDF
        den === nothing && error("PDF must be first precomputed.")
        itp_u = interpolate(den.values, BSpline(Quadratic(Line(OnGrid()))))
        itp = scale(itp_u, den.xs)
        etp = extrapolate(itp, 0.)
        return etp.itp(x)
    elseif method == :naive
        return pdf(ke, x)
    else
        error("Method not available.")
    end
end

# custom broadcast prepares for interpolation only once for all xs
function Base.Broadcast.broadcasted(::typeof(pdf), ke::UnivariateKernelEstimate, xs, method::Symbol)
        if method == :precomputed
            den = ke.precomputedPDF
            den === nothing && error("PDF must be first precomputed.")
            itp_u = interpolate(den.values, BSpline(Quadratic(Line(OnGrid()))))
            itp = scale(itp_u, den.xs)
            etp = extrapolate(itp, 0.)
        return etp.itp.(xs)
    elseif method == :naive
        return pdf.(ke, x)
    else
        error("Method not available.")
    end
end


# it is possible to add cdf(ke, :precomputed)

# it is possibble to add cf(ke, t)