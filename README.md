# KernelDensity.jl

[![CI](https://github.com/JuliaStats/KernelDensity.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaStats/KernelDensity.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/github/JuliaStats/KernelDensity.jl/graph/badge.svg?token=Pvge67IhU8)](https://codecov.io/github/JuliaStats/KernelDensity.jl)
[![Stable Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliastats.org/KernelDensity.jl/stable/)

Kernel density estimators for Julia.

## Usage

### Univariate
The main accessor function is `kde`:

```
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
* `weights`: A vector of weights for each observation. Can be one of `UniformWeights` (the default, from this package), an `AbstractVector` of real numbers, or a `StatsBase.Weights` vector.

The `UnivariateKDE` object `U` contains gridded coordinates (`U.x`) and the density
estimate (`U.density`). These are typically sufficient for plotting.
A related function

``` kde_lscv(data) ```

will construct a `UnivariateKDE` object, with the bandwidth selected by
least-squares cross validation. It accepts the above keyword arguments, except
`bandwidth`.


There are also some slightly more advanced interfaces:
```
kde(data, midpoints::R) where R<:AbstractRange
```
allows specifying the internal grid to use. Optional keyword arguments are
`kernel` and `bandwidth`.

```
kde(data, dist::Distribution)
```
allows specifying the exact distribution to use as the kernel. Optional
keyword arguments are `boundary` and `npoints`.

```
kde(data, midpoints::R, dist::Distribution) where R<:AbstractRange
```
allows specifying both the distribution and grid.

### Bivariate

The usage mirrors that of the univariate case, except that `data` is now
either a tuple of vectors
```
B = kde((xdata, ydata))
```
or a matrix with two columns
```
B = kde(datamatrix)
```
Similarly, the optional arguments all now take tuple arguments:
e.g. `boundary` now takes a tuple of tuples `((xlo,xhi),(ylo,yhi))`.

The `BivariateKDE` object `B` contains gridded coordinates (`B.x` and `B.y`) and the bivariate density
estimate (`B.density`).

### Interpolation
The KDE objects are stored as gridded density values, with attached
coordinates. These are typically sufficient for plotting (see above), but
intermediate values can be evaluated with the `pdf` method (extended from
Distributions.jl), which by default interpolates through the
[FastInterpolations.jl](https://github.com/ProjectTorreyPines/FastInterpolations.jl) backend.

```julia
pdf(k::UnivariateKDE, x)
pdf(k::BivariateKDE, x, y)
```

where `x` and `y` are real numbers or arrays.

The default is a quadratic, C¹-continuous interpolation; you can select a different
scheme — and its boundary condition — with the `method` keyword:

```julia
import FastInterpolations as FI
pdf(k, x; method = FI.LinearInterp())
pdf(k, x; method = FI.CubicInterp())
pdf(k, x; method = FI.CubicInterp(bc = FI.ZeroCurvBC())) # with a boundary condition of your choice
```

See the [FastInterpolations.jl docs](https://projecttorreypines.github.io/FastInterpolations.jl/stable/boundary-conditions/overview/) for the available boundary conditions.

For repeated calls — or to use the Interpolations.jl backend instead — it is
more efficient to construct an interpolation object once and reuse it:

```julia
ik = InterpKDE(k)       # Interpolations.jl backend
ik = FastInterpKDE(k)   # FastInterpolations.jl backend
pdf(ik, x)
```

- `InterpKDE` ([Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl)) passes any extra arguments to `interpolate`.
- `FastInterpKDE` ([FastInterpolations.jl](https://github.com/ProjectTorreyPines/FastInterpolations.jl)) takes a `method` keyword to select the interpolation scheme.