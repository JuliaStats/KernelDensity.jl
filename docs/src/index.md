```@meta
CurrentModule = KernelDensity
```

# KernelDensity.jl

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
* `weights`: A vector of weights for each observation. Can be one of `UniformWeights` (the default, from this package), an `AbstractVector` of real numbers, or a `StatsBase.Weights` vector.

The `UnivariateKDE` object `U` contains gridded coordinates (`U.x`) and the density
estimate (`U.density`). These are typically sufficient for plotting.
A related function

``` kde_lscv(data) ```

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

### Bivariate

The usage mirrors that of the univariate case, except that `data` is now
either a tuple of vectors
```julia
B = kde((xdata, ydata))
```
or a matrix with two columns
```julia
B = kde(datamatrix)
```
Similarly, the optional arguments all now take tuple arguments:
e.g. `boundary` now takes a tuple of tuples `((xlo,xhi),(ylo,yhi))`.

The `BivariateKDE` object `B` contains gridded coordinates (`B.x` and `B.y`) and the bivariate density
estimate (`B.density`).

#### Marginal densities

`B.density[i, j]` estimates the joint density at `(B.x[i], B.y[j])`.
A row holds x fixed and a column holds y fixed. For example, on a rectangular grid:

```@example marginals
using KernelDensity

x = [-1.0, -0.5, 0.0, 0.5, 1.0]
y = [8.0, 12.0, 9.0, 15.0, 11.0]
B = kde((x, y); npoints=(128, 64), boundary=((-4.0, 4.0), (-5.0, 30.0)))
dx, dy = step(B.x), step(B.y)
fx = vec(sum(B.density; dims=2)) .* dy  # integrate over y at each x
fy = vec(sum(B.density; dims=1)) .* dx  # integrate over x at each y
(length(fx), length(fy))
```

For example, `fy[j]` is a density at `B.y[j]`. To estimate the probability
between two adjacent y coordinates, integrate once more: the trapezoid
approximation is `(fy[j] + fy[j+1]) * dy / 2`.
The same calculation with `fx` and `dx` gives the probability for an x interval.

### Interpolation

The KDE objects are stored as gridded density values, with attached
coordinates. These are typically sufficient for plotting (see above), but
intermediate values can be interpolated using the
[Interpolations.jl](https://github.com/tlycken/Interpolations.jl) package via the `pdf` method
(extended from Distributions.jl).

```julia
pdf(k::UnivariateKDE, x)
pdf(k::BivariateKDE, x, y)
```

where `x` and `y` are real numbers or arrays.

If you are making multiple calls to `pdf`, it will be more efficient to
construct an intermediate `InterpKDE` to store the interpolation structure:

```julia
ik = InterpKDE(k)
pdf(ik, x)
```

`InterpKDE` will pass any extra arguments to `interpolate`.

## API Reference

```@autodocs
Modules = [KernelDensity]
```