# KDE.jl

[![Build Status](https://travis-ci.org/simonbyrne/KDE.jl.png)](https://travis-ci.org/simonbyrne/KDE.jl)

Kernel density estimators for julia.

## Usage

### Univariate
The main accessor function is `kde`:
```
kde(data)
```
will construct a `UnivariateKDE` object from the data. The optional keyword arguments are
* `endpoints`: the lower and upper limits of the kde as a tuple. Due to the
  fourier transforms used internally, there should be sufficient spacing to
  prevent wrap-around at the boundaries.
* `npoints`: the number of interpolation points to use. The function uses
  fast fourier transforms (FFTs) internally, so for optimal efficiency this
  should be a power of 2 (default = 2048).
* `kernel`: the distributional family to use as the kernel (default =
  `Normal`). To add your own kernel, extend the internal `kernel_dist` function.
* `bandwidth`: the bandwidth of the kernel. Default is to use Silverman's
  rule.

There are also some slightly more advanced interfaces:
```
kde(data, midpoints::Range)
```
allows specifying the internal grid to use. Optional keyword arguments are
`kernel` and `bandwidth`.

```
kde(data, dist::Distribution)
```
allows specifying the exact distribution to use as the kernel. Optional
keyword arguments are `endpoints` and `npoints`.

```
kde(data, midpoints::Range, dist::Distribution)
```
allows specifying both the distribution and grid.

## To do

* Use an in-place FFT.
* Spline interpolation
* Bias correction
* Improve bandwidth selection
