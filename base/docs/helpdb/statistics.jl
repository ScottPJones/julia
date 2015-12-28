# This file is a part of Julia. License is MIT: http://julialang.org/license

# statistics

"""
    mean(v[, region])

Compute the mean of whole array `v`, or optionally along the dimensions in `region`. Note:
Julia does not ignore `NaN` values in the computation. For applications requiring the
handling of missing data, the `DataArray` package is recommended.
"""
mean

"""
    var(v[, region])

Compute the sample variance of a vector or array `v`, optionally along dimensions in
`region`. The algorithm will return an estimator of the generative distribution's variance
under the assumption that each entry of `v` is an IID drawn from that generative
distribution. This computation is equivalent to calculating `sumabs2(v - mean(v)) /
(length(v) - 1)`. Note: Julia does not ignore `NaN` values in the computation. For
applications requiring the handling of missing data, the `DataArray` package is recommended.
"""
var

"""
    varm(v, m)

Compute the sample variance of a vector `v` with known mean `m`. Note: Julia does not ignore
`NaN` values in the computation.
"""
varm

"""
    stdm(v, m)

Compute the sample standard deviation of a vector `v` with known mean `m`. Note: Julia does
not ignore `NaN` values in the computation.
"""
stdm

"""
    std(v[, region])

Compute the sample standard deviation of a vector or array `v`, optionally along dimensions
in `region`. The algorithm returns an estimator of the generative distribution's standard
deviation under the assumption that each entry of `v` is an IID drawn from that generative
distribution. This computation is equivalent to calculating `sqrt(sum((v - mean(v)).^2) /
(length(v) - 1))`. Note: Julia does not ignore `NaN` values in the computation. For
applications requiring the handling of missing data, the `DataArray` package is recommended.
"""
std

"""
    histrange(v, n)

Compute *nice* bin ranges for the edges of a histogram of `v`, using approximately `n` bins.
The resulting step sizes will be 1, 2 or 5 multiplied by a power of 10. Note: Julia does not
ignore `NaN` values in the computation.
"""
histrange

"""
    hist!(counts, v, e) -> e, counts

Compute the histogram of `v`, using a vector/range `e` as the edges for the bins. This
function writes the resultant counts to a pre-allocated array `counts`.
"""
hist!

"""
    hist(v, e) -> e, counts

Compute the histogram of `v` using a vector/range `e` as the edges for the bins. The result
will be a vector of length `length(e) - 1`, such that the element at location `i` satisfies
`sum(e[i] .< v .<= e[i+1])`. Note: Julia does not ignore `NaN` values in the computation.
"""
hist(v,e)

"""
    hist(v, n) -> e, counts

Compute the histogram of `v`, optionally using approximately `n` bins. The return values are
a range `e`, which correspond to the edges of the bins, and `counts` containing the number
of elements of `v` in each bin. Note: Julia does not ignore `NaN` values in the computation.
"""
hist(v,n::Int=?)

"""
    hist2d!(counts, M, e1, e2) -> (e1, e2, counts)

Compute a "2d histogram" with respect to the bins delimited by the edges given in `e1` and
`e2`. This function writes the results to a pre-allocated array `counts`.
"""
hist2d!

"""
    hist2d(M, e1, e2) -> (edge1, edge2, counts)

Compute a "2d histogram" of a set of N points specified by N-by-2 matrix `M`. Arguments `e1`
and `e2` are bins for each dimension, specified either as integer bin counts or vectors of
bin edges. The result is a tuple of `edge1` (the bin edges used in the first dimension),
`edge2` (the bin edges used in the second dimension), and `counts`, a histogram matrix of
size `(length(edge1)-1, length(edge2)-1)`. Note: Julia does not ignore `NaN` values in the
computation.
"""
hist2d
