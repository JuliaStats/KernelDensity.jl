# KernelDensity.jl

[![CI](https://github.com/JuliaStats/KernelDensity.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaStats/KernelDensity.jl/actions/workflows/CI.yml)
[![Coverage Status](https://coveralls.io/repos/github/JuliaStats/KernelDensity.jl/badge.svg)](https://coveralls.io/github/JuliaStats/KernelDensity.jl)

Kernel density estimators for Julia.

## Usage

### Univariate
The main accessor function is `kde`:

```julia
U = kde(data)
```

will construct a `UnivariateKDE` object from the real vector `data`. The
optional keyword arguments are
* `boundary`: the lower and upper limits of the kde as a tuple. Due to the
  fourier transforms used internally, there should be sufficient spacing to
  prevent wrap-around at the boundaries.
* `npoints`: the number of interpolation points to use. The function uses
  fast Fourier transforms (FFTs) internally, so for optimal efficiency this
  should be a power of 2 (default = 2048).
* `kernel`: the distributional family from
  [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) to use as
  the kernel (default = `Normal`). To add your own kernel, extend the internal
  `kernel_dist` function.
* `bandwidth`: the bandwidth of the kernel. Default is to use Silverman's
  rule.

The `UnivariateKDE` object `U` contains gridded coordinates (`U.x`) and the density
estimate (`U.density`). These are typically sufficient for plotting.
A related function

```julia
kde_lscv(data)
```

will construct a `UnivariateKDE` object, with the bandwidth selected by
least-squares cross validation. It accepts the above keyword arguments, except
`bandwidth`.


There are also some slightly more advanced interfaces:
```julia
kde(data, midpoints::R) where R<:AbstractRange
```
allows specifying the internal grid to use. Optional keyword arguments are
`kernel` and `bandwidth`.

```julia
kde(data, dist::Distribution)
```
allows specifying the exact distribution to use as the kernel. Optional
keyword arguments are `boundary` and `npoints`.

```julia
kde(data, midpoints::R, dist::Distribution) where R<:AbstractRange
```
allows specifying both the distribution and grid.

### Multivariate

The usage mirrors that of the univariate case, except that `data` is now
either an `N`-tuple of vectors
```julia
M = kde((xdata, ydata, zdata))
```
or a matrix with `N` columns
```julia
M = kde(datamatrix)
```
Similarly, the optional arguments all now take `N`-tuple arguments:
e.g. `boundary` now takes a `N`-tuple of tuples `((xlo, xhi), (ylo, yhi), (zlo, zhi))`.

The `MultivariateKDE` object `M` contains gridded coordinates (`M.ranges`, an
`N`-tuple of `AbstractRange`s) and the multivariate density estimate (`M.density`).

### Bi- and Trivariate

Special type definitions exist for the bi- and trivariate case:
```julia
const BivariateKDE = MultivariateKDE{2, ...}
const TrivariateKDE = MultivariateKDE{3, ...}
```

Their contained gridded coordinates can be accessed as `B.x` and `B.y` for
`B isa BivariateKDE` and `T.x`, `T.y`, and `T.z` for `T isa TrivariateKDE`,
respectively.

### Interpolation

The KDE objects are stored as gridded density values, with attached
coordinates. These are typically sufficient for plotting (see above), but
intermediate values can be interpolated using the
[Interpolations.jl](https://github.com/tlycken/Interpolations.jl) package via the `pdf` method
(extended from Distributions.jl).

```julia
pdf(k::UnivariateKDE, x)
pdf(k::BivariateKDE, x, y)
pdf(k::TrivariateKDE, x, y, z)
pdf(k::MultivariateKDE, x...)
```

where `x`, `y`, and `z` are real numbers or arrays.

If you are making multiple calls to `pdf`, it will be more efficient to
construct an intermediate `InterpKDE` to store the interpolation structure:

```julia
ik = InterpKDE(k)
pdf(ik, x)
```

`InterpKDE` will pass any extra arguments to `interpolate`.
