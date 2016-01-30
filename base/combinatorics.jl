# This file is a part of Julia. License is MIT: http://julialang.org/license

# Factorials

const _fact_table64 =
    Int64[1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600,6227020800,
          87178291200,1307674368000,20922789888000,355687428096000,6402373705728000,
          121645100408832000,2432902008176640000]

# Break up parsing of UInt128 constants to avoid using BigInt
@inline pack128(hi::UInt64, lo::UInt64) = UInt128(hi)<<64 | lo

const _fact_table128 =
    UInt128[0x00000000_00000001, 0x00000000_00000002, 0x00000000_00000006, 0x00000000_00000018,
            0x00000000_00000078, 0x00000000_000002d0, 0x00000000_000013b0, 0x00000000_00009d80,
            0x00000000_00058980, 0x00000000_00375f00, 0x00000000_02611500, 0x00000000_1c8cfc00,
            0x00000001_7328cc00, 0x00000014_4c3b2800, 0x00000130_77775800, 0x00001307_77758000,
            0x0001437e_eecd8000, 0x0016beec_ca730000, 0x01b02b93_06890000, 0x21c3677c_82b40000,
	    pack128(0x00000000_00000002, 0xc5077d36_b8c40000),
	    pack128(0x00000000_0000003c, 0xeea4c2b3_e0d80000),
            pack128(0x00000000_00000579, 0x70cd7e29_33680000),
	    pack128(0x00000000_00008362, 0x9343d3dc_d1c00000),
            pack128(0x00000000_000cd4a0, 0x619fb090_7bc00000),
	    pack128(0x00000000_014d9849, 0xea37eeac_91800000),
            pack128(0x00000000_232f0fcb, 0xb3e62c33_58800000),
	    pack128(0x00000003_d925ba47, 0xad2cd59d_ae000000),
            pack128(0x0000006f_99461a1e, 0x9e1432dc_b6000000),
	    pack128(0x00000d13_f6370f96, 0x865df5dd_54000000),
            pack128(0x0001956a_d0aae33a, 0x4560c5cd_2c000000),
	    pack128(0x0032ad5a_155c6748, 0xac18b9a5_80000000),
            pack128(0x0688589c_c0e9505e, 0x2f2fee55_80000000),
	    pack128(0xde1bc4d1_9efcac82, 0x445da75b_00000000)]

function factorial_lookup(n::Integer, table, lim)
    n < 0 && throw(DomainError())
    n > lim && throw(OverflowError())
    n == 0 && return one(n)
    @inbounds f = table[n]
    return oftype(n, f)
end

factorial(n::Int128) = factorial_lookup(n, _fact_table128, 33)
factorial(n::UInt128) = factorial_lookup(n, _fact_table128, 34)
factorial(n::Union{Int64,UInt64}) = factorial_lookup(n, _fact_table64, 20)

if Int === Int32
    factorial(n::Union{Int8,UInt8,Int16,UInt16}) = factorial(Int32(n))
    factorial(n::Union{Int32,UInt32}) = factorial_lookup(n, _fact_table64, 12)
else
    factorial(n::Union{Int8,UInt8,Int16,UInt16,Int32,UInt32}) = factorial(Int64(n))
end

function gamma(n::Union{Int8,UInt8,Int16,UInt16,Int32,UInt32,Int64,UInt64})
    n < 0 && throw(DomainError())
    n == 0 && return Inf
    n <= 2 && return 1.0
    n > 20 && return gamma(Float64(n))
    @inbounds return Float64(_fact_table64[n-1])
end


# Basic functions for working with permutations

"""
    isperm(v) -> Bool

Returns `true` if `v` is a valid permutation.
"""
function isperm(A)
    n = length(A)
    used = falses(n)
    for a in A
        (0 < a <= n) && (used[a] $= true) || return false
    end
    true
end

function permute!!{T<:Integer}(a, p::AbstractVector{T})
    count = 0
    start = 0
    while count < length(a)
        ptr = start = findnext(p, start+1)
        temp = a[start]
        next = p[start]
        count += 1
        while next != start
            a[ptr] = a[next]
            p[ptr] = 0
            ptr = next
            next = p[next]
            count += 1
        end
        a[ptr] = temp
        p[ptr] = 0
    end
    a
end

"""
    permute!(v, p)

Permute vector `v` in-place, according to permutation `p`. No checking is done
to verify that `p` is a permutation.

To return a new permutation, use `v[p]`. Note that this is generally faster than
`permute!(v,p)` for large vectors.
"""
permute!(a, p::AbstractVector) = permute!!(a, copymutable(p))

function ipermute!!{T<:Integer}(a, p::AbstractVector{T})
    count = 0
    start = 0
    while count < length(a)
        start = findnext(p, start+1)
        temp = a[start]
        next = p[start]
        count += 1
        while next != start
            temp_next = a[next]
            a[next] = temp
            temp = temp_next
            ptr = p[next]
            p[next] = 0
            next = ptr
            count += 1
        end
        a[next] = temp
        p[next] = 0
    end
    a
end

"""
    ipermute!(v, p)

Like `permute!`, but the inverse of the given permutation is applied.
"""
ipermute!(a, p::AbstractVector) = ipermute!!(a, copymutable(p))

"""
    invperm(v)

Return the inverse permutation of `v`
"""
function invperm(a::AbstractVector)
    b = zero(a) # similar vector of zeros
    n = length(a)
    @inbounds for (i, j) in enumerate(a)
        ((1 <= j <= n) && b[j] == 0) ||
            throw(ArgumentError("argument is not a permutation"))
        b[j] = i
    end
    b
end

#XXX This function should be moved to Combinatorics.jl but is currently used by Base.DSP.
"""
    nextprod([k_1,k_2,...], n)

Next integer not less than `n` that can be written as ``\\prod k_i^{p_i}`` for integers
``p_1``, ``p_2``, etc.
"""
function nextprod(a::Vector{Int}, x)
    if x > typemax(Int)
        throw(ArgumentError("unsafe for x > typemax(Int), got $x"))
    end
    k = length(a)
    v = ones(Int, k)                  # current value of each counter
    mx = [nextpow(ai,x) for ai in a]  # maximum value of each counter
    v[1] = mx[1]                      # start at first case that is >= x
    p::widen(Int) = mx[1]             # initial value of product in this case
    best = p
    icarry = 1

    while v[end] < mx[end]
        if p >= x
            best = p < best ? p : best  # keep the best found yet
            carrytest = true
            while carrytest
                p = div(p, v[icarry])
                v[icarry] = 1
                icarry += 1
                p *= a[icarry]
                v[icarry] *= a[icarry]
                carrytest = v[icarry] > mx[icarry] && icarry < k
            end
            if p < x
                icarry = 1
            end
        else
            while p < x
                p *= a[1]
                v[1] *= a[1]
            end
        end
    end
    # might overflow, but want predictable return type
    return mx[end] < best ? Int(mx[end]) : Int(best)
end
