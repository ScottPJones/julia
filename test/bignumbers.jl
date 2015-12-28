# This file is a part of Julia. License is MIT: http://julialang.org/license

if Base.BUILD_BIGINT # SPJ!!! Could some of these be rewritten not to use BigInt?
@test muladd(big(1),2,3) == big(1)*2+3
@test muladd(big(1//1),2,3) == big(1//1)*2+3
@test muladd(big(1.0),2,3) == big(1.0)*2+3

let x = BigInt(7)^77
    @test fma(x-1, x-2, x-3) == (x-1) * (x-2) + (x-3)
    @test (fma((x-1)//(x-2), (x-3)//(x-4), (x-5)//(x-6)) ==
           (x-1)//(x-2) * (x-3)//(x-4) + (x-5)//(x-6))
end

let eps = 1//BigInt(2)^30, one_eps = 1+eps,
    eps64 = Float64(eps), one_eps64 = Float64(one_eps)
    @test eps64 == Float64(eps)
    @test rationalize(BigInt, eps64, tol=0) == eps
    @test one_eps64 == Float64(one_eps)
    @test rationalize(BigInt, one_eps64, tol=0) == one_eps
    @test one_eps64 * one_eps64 - 1 != Float64(one_eps * one_eps - 1)
    @test fma(one_eps64, one_eps64, -1) == Float64(one_eps * one_eps - 1)
end

let eps = 1//BigInt(2)^15, one_eps = 1+eps,
    eps32 = Float32(eps), one_eps32 = Float32(one_eps)
    @test eps32 == Float32(eps)
    @test rationalize(BigInt, eps32, tol=0) == eps
    @test one_eps32 == Float32(one_eps)
    @test rationalize(BigInt, one_eps32, tol=0) == one_eps
    @test one_eps32 * one_eps32 - 1 != Float32(one_eps * one_eps - 1)
    @test fma(one_eps32, one_eps32, -1) == Float32(one_eps * one_eps - 1)
end

let eps = 1//BigInt(2)^7, one_eps = 1+eps,
    eps16 = Float16(Float32(eps)), one_eps16 = Float16(Float32(one_eps))
    @test eps16 == Float16(Float32(eps))
    # Currently broken in Julia -- enable when "rationalize" is fixed;
    # see <https://github.com/JuliaLang/julia/issues/9897>
    # @test rationalize(BigInt, eps16, tol=0) == eps
    @test one_eps16 == Float16(Float32(one_eps))
    # @test rationalize(BigInt, one_eps16, tol=0) == one_eps
    @test one_eps16 * one_eps16 - 1 != Float16(Float32(one_eps * one_eps - 1))
    @test (fma(one_eps16, one_eps16, -1) ==
           Float16(Float32(one_eps * one_eps - 1)))
end

let eps = 1//BigInt(2)^200, one_eps = 1+eps,
    eps256 = BigFloat(eps), one_eps256 = BigFloat(one_eps)
    @test eps256 == BigFloat(eps)
    @test rationalize(BigInt, eps256, tol=0) == eps
    @test one_eps256 == BigFloat(one_eps)
    @test rationalize(BigInt, one_eps256, tol=0) == one_eps
    @test one_eps256 * one_eps256 - 1 != BigFloat(one_eps * one_eps - 1)
    @test fma(one_eps256, one_eps256, -1) == BigFloat(one_eps * one_eps - 1)
end

# muladd

let eps = 1//BigInt(2)^30, one_eps = 1+eps,
    eps64 = Float64(eps), one_eps64 = Float64(one_eps)
    @test eps64 == Float64(eps)
    @test one_eps64 == Float64(one_eps)
    @test one_eps64 * one_eps64 - 1 != Float64(one_eps * one_eps - 1)
    @test isapprox(muladd(one_eps64, one_eps64, -1),
                   Float64(one_eps * one_eps - 1))
end

let eps = 1//BigInt(2)^15, one_eps = 1+eps,
    eps32 = Float32(eps), one_eps32 = Float32(one_eps)
    @test eps32 == Float32(eps)
    @test one_eps32 == Float32(one_eps)
    @test one_eps32 * one_eps32 - 1 != Float32(one_eps * one_eps - 1)
    @test isapprox(muladd(one_eps32, one_eps32, -1),
                   Float32(one_eps * one_eps - 1))
end

let eps = 1//BigInt(2)^7, one_eps = 1+eps,
    eps16 = Float16(Float32(eps)), one_eps16 = Float16(Float32(one_eps))
    @test eps16 == Float16(Float32(eps))
    @test one_eps16 == Float16(Float32(one_eps))
    @test one_eps16 * one_eps16 - 1 != Float16(Float32(one_eps * one_eps - 1))
    @test isapprox(muladd(one_eps16, one_eps16, -1),
                   Float16(Float32(one_eps * one_eps - 1)))
end

# large integer literals
@test isa(170141183460469231731687303715884105728,BigInt)
@test isa(0170141183460469231731687303715884105728,BigInt)
@test isa(-170141183460469231731687303715884105729,BigInt)

# exponentiating with a negative base
@test -170141183460469231731687303715884105728^2 ==
    -(170141183460469231731687303715884105728^2)

# numeric literal coefficients
let x = 10
    @test 170141183460469231731687303715884105728x ==
        1701411834604692317316873037158841057280
end

@test 170141183460469231731687303715884105728(10) ==
    1701411834604692317316873037158841057280

# verify type stability with integer (x is negative)
@test eltype(copysign(-1,BigInt(1))) <: Integer
@test eltype(copysign(-BigInt(1),1)) <: Integer
@test eltype(copysign(-BigInt(1),1.0)) <: Integer
@test eltype(copysign(-BigInt(1),1//2)) <: Integer
@test eltype(copysign(-BigInt(1),BigInt(1))) <: Integer
@test eltype(copysign(-1,-BigInt(1))) <: Integer
@test eltype(copysign(-BigInt(1),-1)) <: Integer
@test eltype(copysign(-BigInt(1),-1.0)) <: Integer
@test eltype(copysign(-BigInt(1),-1//2)) <: Integer
@test eltype(copysign(-BigInt(1),-BigInt(1))) <: Integer

# verify type stability with integer (x is positive)
@test eltype(copysign(1,BigInt(1))) <: Integer
@test eltype(copysign(BigInt(1),1)) <: Integer
@test eltype(copysign(BigInt(1),1.0)) <: Integer
@test eltype(copysign(BigInt(1),1//2)) <: Integer
@test eltype(copysign(BigInt(1),BigInt(1))) <: Integer
@test eltype(copysign(1,-BigInt(1))) <: Integer
@test eltype(copysign(BigInt(1),-1)) <: Integer
@test eltype(copysign(BigInt(1),-1.0)) <: Integer
@test eltype(copysign(BigInt(1),-1//2)) <: Integer
@test eltype(copysign(BigInt(1),-BigInt(1))) <: Integer

# verify type stability with real (x is negative)
@test eltype(copysign(-1.0,BigInt(1))) <: Real
@test eltype(copysign(-1.0,-BigInt(1))) <: Real

# Verify type stability with real (x is positive)
@test eltype(copysign(1.0,BigInt(1))) <: Real
@test eltype(copysign(1.0,-BigInt(1))) <: Real

# Verify type stability with rational (x is negative)
@test eltype(copysign(-1//2,BigInt(1))) <: Rational
@test eltype(copysign(-1//2,-BigInt(1))) <: Rational

# Verify type stability with rational (x is positive)
@test eltype(copysign(-1//2,BigInt(1))) <: Rational
@test eltype(copysign(-1//2,-BigInt(1))) <: Rational

@test realmin() != 1//(BigInt(2)^1022+1)
@test realmin() == 1//(BigInt(2)^1022)
@test realmin() != 1//(BigInt(2)^1022-1)
@test realmin()/2 != 1//(BigInt(2)^1023+1)
@test realmin()/2 == 1//(BigInt(2)^1023)
@test realmin()/2 != 1//(BigInt(2)^1023-1)
@test nextfloat(0.0) != 1//(BigInt(2)^1074+1)
@test nextfloat(0.0) == 1//(BigInt(2)^1074)
@test nextfloat(0.0) != 1//(BigInt(2)^1074-1)

# issue 6712
@test convert(Rational{BigInt},Float64(pi)) == Float64(pi)
@test convert(Rational{BigInt},big(pi)) == big(pi)

# issue 3412
@test convert(Rational{BigInt},0.0) == 0
@test convert(Rational{BigInt},-0.0) == 0
@test convert(Rational{BigInt},zero(BigFloat)) == 0
@test convert(Rational{BigInt},-zero(BigFloat)) == 0
@test convert(Rational{BigInt},5e-324) == 5e-324
@test convert(Rational{BigInt},realmin(Float64)) == realmin(Float64)
@test convert(Rational{BigInt},realmax(Float64)) == realmax(Float64)

@test isa(convert(Float64, big(1)//2), Float64)

@test rationalize(BigInt,nextfloat(0.1)) == 300239975158034//3002399751580339
@test rationalize(BigInt,nextfloat(0.1),tol=0.5eps(0.1)) == 379250494936463//3792504949364629
@test rationalize(BigInt,nextfloat(0.1),tol=1.5eps(0.1)) == 1//10
@test rationalize(BigInt,nextfloat(0.1),tol=0) == 7205759403792795//72057594037927936
@test rationalize(BigInt,prevfloat(0.1)) == 1//10
@test rationalize(BigInt,prevfloat(0.1),tol=0) == 7205759403792793//72057594037927936

# factorization of factors > 2^16
@test factor((big(2)^31-1)^2) == Dict(big(2^31-1) => 2)
@test factor((big(2)^31-1)*(big(2)^17-1)) == Dict(big(2^31-1) => 1, big(2^17-1) => 1)

# check gcd and related functions against GMP
for T in (Int32,Int64), ii = -20:20, jj = -20:20
    i::T, j::T = ii, jj
    local d = gcd(i,j)
    @test d >= 0
    @test lcm(i,j) >= 0
    local ib = big(i)
    local jb = big(j)
    @test d == gcd(ib,jb)
    @test lcm(i,j) == lcm(ib,jb)
    @test gcdx(i,j) == gcdx(ib,jb)
    if j == 0
        @test_throws ErrorException invmod(i,j)
        @test_throws ErrorException invmod(ib,jb)
    elseif d == 1
        n = invmod(i,j)
        @test n == invmod(ib,jb)
        @test mod(n*i,j) == mod(1,j)
    end
end

for T in (Int8,Int16,Int32,Int64,Int128,UInt8,UInt16,UInt32,UInt64,UInt128)
    @test_throws InexactError T(big(typemax(T))+1)
    @test_throws InexactError T(big(typemin(T))-1)
end

@test digits(BigInt(2)^128, 2) == [zeros(128); 1]

@test_throws InexactError convert(UInt8, big(300))
@test_throws InexactError convert(Int, big(2)^100)
@test_throws InexactError convert(Int16, big(2)^100)

# issue #9611
@test factor(Int128(2)^101+1) == Dict(3=>1,845100400152152934331135470251=>1)
# test second branch, after all small primes in list have been searched
@test factor(10009 * Int128(1000000000000037)) == Dict(10009=>1,1000000000000037=>1)

# Issue #9618: errors thrown by large exponentiations
@test_throws DomainError big(2)^-(big(typemax(UInt))+1)
@test_throws OverflowError big(2)^(big(typemax(UInt))+1)
@test 0==big(0)^(big(typemax(UInt))+1)

# issue #10311
let n = 1
    @test n//n + n//big(n)*im == 1//1 + 1//1*im
end

# BigInt - (small negative) is tricky because gmp only has gmpz_sub_ui
@test big(-200) - Int8(-128) == -72

# check powermod function against few types (in particular [U]Int128 and BigInt)

# with m==1 should give 0
@test powermod(1,0,big(1)) == 0
@test powermod(1,0,big(-1)) == 0
# divide by zero error
@test_throws DivideError powermod(1,0,big(0))
# negative power domain error
@test_throws DomainError powermod(1,-2,big(1))

for i = -100:100
    @test nextpow2(i) == nextpow2(big(i))
    @test prevpow2(i) == prevpow2(big(i))
end

@test widen(BigInt) === BigInt
@test typeof(widemul(Int128(1),UInt128(1))) == BigInt
@test typeof(widemul(UInt128(1),Int128(1))) == BigInt
end

if Base.BUILD_BIGFLT
# issue 3412
@test convert(Rational,zero(BigFloat)) == 0
@test convert(Rational,-zero(BigFloat)) == 0

# issue #6365
for T in (Float32, Float64)
    for i = 9007199254740992:9007199254740996
        @test T(i) == T(BigFloat(i))
        @test T(-i) == T(BigFloat(-i))
        for r in (RoundNearest,RoundUp,RoundDown,RoundToZero)
            @test T(i,r) == T(BigFloat(i),r)
            @test T(-i,r) == T(BigFloat(-i),r)
        end
    end
end

@test prevfloat(big(pi)) < pi
@test nextfloat(big(pi)) > pi
@test !(prevfloat(big(pi)) > pi)
@test !(nextfloat(big(pi)) < pi)

## Note: this should change to e.g. Float128 at some point
@test widen(Float64) === BigFloat

let x = big(-0.0)
    @test signbit(x) && !signbit(abs(x))
end

if Base.BUILD_BIGINT
@test rationalize(BigInt,nextfloat(parse(BigFloat,"0.1")),tol=1.5eps(big(0.1))) == 1//10
@test rationalize(BigInt,prevfloat(parse(BigFloat,"0.1"))) == 1//10
@test rationalize(BigInt,nextfloat(parse(BigFloat,"0.1")),tol=0) == 46316835694926478169428394003475163141307993866256225615783033603165251855975//463168356949264781694283940034751631413079938662562256157830336031652518559744

for (d,B) in ((4//2+1im,Rational{BigInt}),(3.0+1im,BigFloat),(2+1im,BigInt))
    @test typeof(big(d)) == Complex{B}
    @test big(d) == d
    @test typeof(big([d])) == Vector{Complex{B}}
    @test big([d]) == [d]
end
end
end

