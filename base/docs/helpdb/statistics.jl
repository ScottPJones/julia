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
