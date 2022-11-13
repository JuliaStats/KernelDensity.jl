"""
$(TYPEDEF)

Store both grid and density for KDE over the real line.

Reading the fields directly is part of the API, and

```julia
sum(density) * prod(step.(ranges)) ≈ 1
```

# Fields

$(FIELDS)
"""
struct MultivariateKDE{N, R <: NTuple{N, AbstractRange}} <: AbstractKDE
    "Coordinates of gridpoints for evaluating density"
    ranges::R
    "Kernel density at corresponding gridpoints."
    density::Array{Float64, N}
end

function kernel_dist(::Type{D}, w::NTuple{N, Real}) where {N, D <: UnivariateDistribution}
    kernel_dist.(D, w)
end

@generated function kernel_dist(::Type{Ds}, w::NTuple{N, Real}) where {N, Ds <: NTuple{N, UnivariateDistribution}}
    quote
        kernel_dist.(
            # Ds.parameters is of type svec which does not have a compile time
            # known length. By splatting it into a tuple at compile time, proper
            # type inference is possible.
            ($(Ds.parameters...),),
            w
        )
    end
end

const DataTypeOrUnionAll = Union{DataType, UnionAll}

# this function provided for backwards compatibility, though it doesn't have the type restrictions
# to ensure that the given tuple only contains univariate distributions
function kernel_dist(d::NTuple{N, DataTypeOrUnionAll}, w::NTuple{N, Real}) where N
    kernel_dist.(d, w)
end

# TODO: there are probably better choices.
function default_bandwidth(data::NTuple{N, RealVector}) where N
    default_bandwidth.(data)
end

@generated function interpolate_in_hypercube!(
    grid::Array{Float64, N},
    midpoints::NTuple{N, AbstractRange},
    coords::NTuple{N, Real},
    high_idcs::NTuple{N, Int},
    low_idcs::NTuple{N, Int},
    weight::Real
) where N

    interpolate_in_hypercube!_impl(Val(N))

end

function interpolate_in_hypercube!_impl(::Val{N}) where N
    # We need to distribute the `weight` of one data point to all vertices of
    # the hypercube that surrounds it in the grid.
    # These vertices can be described using all N-bit bitstrings.
    # A zero bit then means that we are looking at the lower index in a specific
    # dimension and a one means we look at the higher one.
    side_to_symbol(side) = side == 0 ? :low_idcs : :high_idcs
    updates = quote end
    sides = zeros(Int, N)

    for vertex_number in 0 : 2^N - 1
        # Get the bitstring representation of the vertex.
        digits!(sides, vertex_number, base = 2)
        # Reverse it such that we will iterate the last dimension the fastest
        # to adhere to Julia's "column"-major memory layout of arrays.
        reverse!(sides)

        indices = [:($(side_to_symbol(side))[$idx]) for (side, idx) in zip(sides, 1:N)]

        # The interpolation coefficients are chosen such that the partial weight
        # on one vertex is proportional to the distance of the data point to the
        # opposite vertex for each dimension respectively.
        factors = [
            :(midpoints[$idx][$(side_to_symbol(1-side))[$idx]] - coords[$idx])
            for (side, idx) in zip(sides, 1:N)
        ]
        # The previous formula is only almost correct.
        # Whenever we calculate the distance to a vertex with a lower index, we
        # introduce a multiplicative error of -1.
        # This is mitigated by the following correction:
        sign_correction = (-1)^sum(sides)

        updates = quote
            $updates

            @inbounds grid[$(indices...)] += *($sign_correction, weight, $(factors...))
        end
    end

    updates
end

# tabulate data for kde
function tabulate(
    data::NTuple{N, RealVector},
    midpoints::NTuple{N, AbstractRange},
    weights::Weights = default_weights(data)
) where N

    ndata = length(first(data))
    all(==(ndata) ∘ length, data) || error("data vectors must be of same length")

    range_lengths = length.(midpoints)
    range_steps = step.(midpoints)

    # Set up a grid for discretized data
    grid = zeros(Float64, range_lengths...)
    ainc = 1.0 / (sum(weights) * prod(range_steps)^2)

    # weighted discretization (cf. Jones and Lotwick)
    for i in 1:ndata
        coords = getindex.(data, i)
        high_idcs = searchsortedfirst.(midpoints, coords)
        low_idcs = high_idcs .- 1
        if all(1 .<= low_idcs .<= range_lengths .- 1)
            interpolate_in_hypercube!(grid, midpoints, coords, high_idcs, low_idcs, weights[i])
        end
    end
    # Normalize for interpolation coefficients and weights
    grid .*= ainc

    # returns an un-convolved KDE
    MultivariateKDE(midpoints, grid)
end

# convolution with product distribution of N univariate distributions
function conv(k::MultivariateKDE{N, R}, dists::NTuple{N, UnivariateDistribution}) where {N, R}
    # Transform to Fourier basis
    ft = rfft(k.density)

    # Convolve fft with characteristic function of kernel
    cs = -twoπ ./ (maximum.(k.ranges) .- minimum.(k.ranges))
    for idx in CartesianIndices(ft)
        pos = Tuple(idx) .- 1
        pos = min.(pos, size(k.density) .- pos)
        ft[idx] *= prod(cf.(dists, pos .* cs))
    end

    # Invert the Fourier transform to get the KDE
    density = irfft(ft, size(k.density, 1))

    density .= max.(0., density)

    MultivariateKDE(k.ranges, density)
end

default_weights(data::NTuple{N, RealVector}) where N = UniformWeights(length(data[1]))

function kde(
    data::NTuple{N, RealVector},
    weights::Weights,
    midpoints::NTuple{N, AbstractRange},
    dist::NTuple{N, UnivariateDistribution}
) where N

    k = tabulate(data, midpoints, weights)
    conv(k, dist)
end

function kde(
    data::NTuple{N, RealVector},
    dist::NTuple{N, UnivariateDistribution};
    boundary::NTuple{N, Tuple{Real, Real}} =
        kde_boundary.(data, std.(dist)),
    npoints::NTuple{N, Int} = ntuple(_ -> 256, Val(N)),
    weights::Weights = default_weights(data)
) where N

    midpoints = kde_range.(boundary, npoints)
    kde(data, weights, midpoints, dist)
end

function kde(
    data::NTuple{N, RealVector},
    midpoints::NTuple{N, AbstractRange};
    bandwidth = default_bandwidth(data),
    kernel = Normal,
    weights::Weights = default_weights(data)
) where N

    dist = kernel_dist(kernel, bandwidth)
    kde(data, weights, midpoints, dist)
end

function kde(
    data::NTuple{N, RealVector};
    bandwidth = default_bandwidth(data),
    kernel = Normal,
    boundary::NTuple{N, Tuple{Real, Real}} = kde_boundary.(data, bandwidth),
    npoints::NTuple{N, Int} = ntuple(_ -> 256, Val(N)),
    weights::Weights = default_weights(data)
) where N

    dist = kernel_dist(kernel, bandwidth)
    midpoints = kde_range.(boundary, npoints)

    kde(data, weights, midpoints, dist)
end

# matrix data
function kde(data::RealMatrix, args...; kwargs...)
    kde(
        ntuple(i -> view(data, :, i), size(data, 2)),
        args...; kwargs...
    )
end
