using Documenter
using KernelDensity

makedocs(
    sitename = "KernelDensity.jl",
    modules = [KernelDensity],
)
deploydocs(
    repo = "github.com/JuliaStats/KernelDensity.jl.git",
    push_preview = true,
)
