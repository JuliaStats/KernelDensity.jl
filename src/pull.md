The current capabilities of `KernelDensity` are quite limited, especially if we compare them to what is avaiable e.g. in SciPy, and this supposed to be is the main source of kernel estimation tools in Julia statistical ecosystem.  The current interface and type system were made for calculating discretised pdf using Fourier transform and not for much else. This is a proposition of completely new type system and interface.

The main principles of this proposition are:
- The new `KernelEstimate` datatype contains full information about the kernel fit, so a statistician can infer all they possibly want from it.
- `KernelEstimate` is a subtype of `MixtureDistribution` so it fits well the rest of `JuliaStats`, can be much more easily used in other tasks and has a lot of useful default methods already implemented in `Distributions`.
- This type system can be used for kernel estimation in any dimension.
- The discretised values of the pdf are distinct from `KernelEstimate`. However, they can be stored in it and if so, getting pdf from interpolation of it becomes avaiable. If it is not present, different methods for getting pdf can be used (e.g. from the naive sum formula, which for small samples is completely reasonable). It can be also added using `precompute!`.

The current approach does not allow for adaptive kde. You can allow for it with the same design, but for the sake of efficiency and code simplicity it might be better to add something like `AdaptiveKernelEstimate` if such a need arises.

I've curretly implemented methods only for 1D. Translating 2D case should be easy, but it might be best to forget about it and implement any dimension from the get go. But this is better to do that after getting some feedback.

Old interface may, or may not stay for backwards compability.

As it stands this proposition does not add any new numerical methods, but even simply linking it to `Distributions` adds a lot of capabilities to this library, which for example Scipy or Sklearn do not have, making it I hope quite more attractive.