# use scale instead of std
# LocationScale ?
# you can use support to detect 
# argument for boundary?
# give any kernel?
# fix cosine distribution

using Plots

using Distributions
using .KernelDensity

##

points = [1. 2 3 4 5]
cat = Categorical([1/5,1/5,1/5,1/5,1/5])

ds = DiscretisedPDF(fill(1/20,20),LinRange(0,1,20))

k = KernelEstimate(points,cat, Normal(0,1),1.0, nothing)
ncomponents(k1)
component(k1,3)
probs(k1)

##

X = randn(10^5)

k1 = kernel_estimate(X, 0.2, Normal)

k2 = kernel_estimate(X, Silverman(), Epanechnikov)

k3 = kernel_estimate(X, LSCV(), Epanechnikov)


precompute!(k2,2048,(-5,5))

pdf(k2, 1)
pdf.(k2, [1, 2, 3])
pdf(k2, 1, :precomputed)
pdf.(k2, [1, 2, 3], :precomputed)

cdf(k2, 1)
mean(k2), var(k2)
quantile(k2, 0.9)


kde = k2.precomputedPDF
plot(kde.xs,kde.values)


