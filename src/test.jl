# use scale instead of std
# LocationScale ?
# you can use support to detect 
# argument for boundary?
# give any kernel?
# fix cosine distribution

# add precompute!

using Distributions
using Plots
using .KernelDensity

#

points = [1. 2 3 4 5]
cat = Categorical([1/5,1/5,1/5,1/5,1/5])

ds = DiscretisedPDF(fill(1/20,20),LinRange(0,1,20))

k = KernelEstimate(points,cat, Normal(0,1),1.0, nothing)

X = randn(10^5)

k1 = kernel_estimate(X,0.2,Normal)

k2 = kernel_estimate(X,Silverman(),Epanechnikov)

precompute!(k2)


kde = k1.precomputedPDF
plot(kde.xs,kde.values)

KernelEstimate(reshape(X,1,length(X)), cat, kernel_dist(Normal,1.0), 1.0)

ncompoments(k)
component(k,3)
probs(k)