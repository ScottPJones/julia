# This file is a part of Julia. License is MIT: http://julialang.org/license

immutable Rational{T<:Integer} <: Real
    num::T
    den::T

    function Rational(num::Integer, den::Integer)
        num == den == zero(T) && throw(ArgumentError("invalid rational: zero($T)//zero($T)"))
        g = den < 0 ? -gcd(den, num) : gcd(den, num)
        new(div(num, g), div(den, g))
    end
end

function divgcd end
function // end
function .// end
function rationalize end
function num end
function den end

immutable Complex{T<:Real} <: Number
    re::T
    im::T
end

typealias Complex32  Complex{Float16}
typealias Complex64  Complex{Float32}
typealias Complex128 Complex{Float64}

function real end
function imag end
function complex end
function isreal end
function isimag end

real(x::Real) = x
imag(x::Real) = zero(x)
real{T<:Real}(::Type{T}) = T
isreal(x::Real) = true
isimag(z::Number) = real(z) == 0

# Add generic function definitions
function log end
function exp end
function sin end
function cos end
function tan end
function sinh end
function cosh end
function tanh end
function asin end
function acos end
function atan end
function asinh end
function acosh end
function atanh end
function sqrt end
function log2 end
function log10 end
function exp2 end
function exp10 end
function expm1 end
function log1p end

function cis end
function reim end
