using Test

tests = [
         "univariate",
         "bivariate",
         "interp",
         ]

@testset "$t" for t in tests
    test_fn = "$t.jl"
    include(test_fn)
end
