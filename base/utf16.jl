# This file is a part of Julia. License is MIT: http://julialang.org/license

#=
@doc """
@brief      Error messages for Unicode / UTF support
""" ->
=#
#=
@enum(UTF_ERR,
      SHORT = 1,
      CONT,
      LONG,
      NOT_LEAD,
      NOT_TRAIL,
      NOT_SURROGATE,
      MISSING_SURROGATE,
      INVALID,
      SURROGATE,
      NULL_TERMINATE,
)
=#

const UTF_ERR_SHORT = 1
const UTF_ERR_CONT  = 2
const UTF_ERR_LONG  = 3
const UTF_ERR_NOT_LEAD = 4
const UTF_ERR_NOT_TRAIL = 5
const UTF_ERR_NOT_SURROGATE = 6
const UTF_ERR_MISSING_SURROGATE = 7
const UTF_ERR_INVALID = 8
const UTF_ERR_SURROGATE = 9
const UTF_ERR_NULL_TERMINATE = 10
const UTF_ERR_MAX = 10

type UnicodeError <: Exception
    msg::AbstractString
end

const errMsgs = [
    "invalid UTF-8 sequence starting at index <<1>> (0x<<2>>) missing one or more continuation bytes)",
    "invalid UTF-8 sequence starting at index <<1>> (0x<<2>> is not a continuation byte)",
    "invalid UTF-8 sequence, overlong encoding starting at index <<1>> (0x<<2>>)",
    "not a leading Unicode surrogate character at index <<1>> (0x<<2>>)",
    "not a trailing Unicode surrogate character at index <<1>> (0x<<2>>)",
    "not a valid Unicode surrogate character at index <<1>> (0x<<2>>",
    "missing trailing Unicode surrogate character after index <<1>> (0x<<2>>)",
    "invalid Unicode character starting at index <<1>> (0x<<2>> > 0x10ffff)",
    "surrogate encoding not allowed in UTF-8 or UTF-32, at index <<1>> (0x<<2>>)",
    "UTF16String data must be NULL-terminated"
]
#=
@doc """
@brief      Throws ArgumentError with information about the specific error, location, and character

@param[in]  errCode::UTF_ERR
@param[in]  errPos:: Integer
@param[in]  errChar::Integer

@throws never returns, always throws ArgumentError
""" ->
=#
function utf_errfunc(errCode::Integer, errPos::Integer, errChar::Integer)
    if errCode < 1 || errCode > UTF_ERR_MAX
        throw(ArgumentError("Invalid error code for Unicode error: $errCode, Pos = $errPos, Char = $errChar"))
    end
    throw(ArgumentError(replace(replace(errMsgs[errCode],"<<1>>",string(errPos)),"<<2>>",hex(errChar))))
    -1
end

#=
@doc """
@brief      Base UTF16String type, has 16-bit NULL termination word after data, native byte order
""" ->
=#
immutable UTF16String <: AbstractString
    data::Vector{UInt16} # includes 16-bit NULL termination after string chars
    function UTF16String(data::Vector{UInt16})
        if length(data) < 1 || data[end] != 0
            utf_errfunc(UTF_ERR_NULL_TERMINATE, 0, 0)
        end
        new(data)
    end
end

#=
@doc """
@brief      Base UTF32String type, has 32-bit NULL termination word after data, native byte order
""" ->
=#
immutable UTF32String <: DirectIndexString
    data::Vector{Char} # includes 32-bit NULL termination after string chars

    function UTF32String(a::Vector{Char})
        if length(a) < 1 || a[end] != Char(0)
            throw(ArgumentError("UTF32String data must be NULL-terminated"))
        end
        new(a)
    end
end

utf16_is_lead(c::UInt16) = (c & 0xfc00) == 0xd800
utf16_is_trail(c::UInt16) = (c & 0xfc00) == 0xdc00
utf16_is_surrogate(c::UInt16) = (c & 0xf800) == 0xd800
utf16_get_supplementary(lead::UInt16, trail::UInt16) = Char(UInt32(lead-0xd7f7)<<10 + trail)

is_surrogate_lead(c::Unsigned) = ((c & 0xfffffc00) == 0xd800)
is_surrogate_trail(c::Unsigned) = ((c & 0xfffffc00) == 0xdc00)
is_surrogate_char(c::Unsigned) = ((c & 0xfffff800) == 0xd800)
is_valid_continuation(c) = ((c & 0xc0) == 0x80)

function length(s::UTF16String)
    d = s.data
    len = length(d) - 1
    len == 0 && return 0
    cnum = 0
    for i = 1:len
        @inbounds cnum += !utf16_is_trail(d[i])
    end
    cnum
end

function endof(s::UTF16String)
    d = s.data
    i = length(d) - 1
    i == 0 && return i
    return utf16_is_surrogate(d[i]) ? i-1 : i
end

function next(s::UTF16String, i::Int)
    ch = s.data[i]
    !utf16_is_surrogate(ch) && return (Char(ch), i+1)
    # check length, account for terminating \0
    i >= (length(s.data)-1) && utf_errfunc(UTF_ERR_MISSING_SURROGATE, i, UInt32(ch))
    !utf16_is_lead(ch) && utf_errfunc(UTF_ERR_NOT_LEAD, i, ch)
    ct = s.data[i+1]
    !utf16_is_trail(ct) && utf_errfunc(UTF_ERR_NOT_TRAIL, i, ch)
    utf16_get_supplementary(ch, ct), i+2
end

function reverseind(s::UTF16String, i::Integer)
    j = length(s.data) - i
    return Base.utf16_is_trail(s.data[j]) ? j-1 : j
end

lastidx(s::UTF16String) = length(s.data) - 1 # s.data includes NULL terminator

function reverse(s::UTF16String)
    d = s.data
    out = similar(d)
    out[end] = 0 # NULL termination
    n = length(d)
    @inbounds for i = 1:n-1
        ch = d[n-i]
        if Base.utf16_is_lead(ch)
            out[i],out[i-1] = out[i-1],ch
        else
            out[i] = ch
        end
    end
    UTF16String(out)
end

#=
function convert(::Type{UTF16String}, str::AbstractString)
    len, cnt2::UInt, cnt3::UInt, cnt4::UInt, flags::UInt = check_string(UTF32String, str)
    buf = Vector{UInt16}(len+cnt4+1)
    i = 0
    @inbounds for ch in str
        c = reinterpret(UInt32, ch)
        if c < 0x10000
            buf[i += 1] = UInt16(c)
        elseif c <= 0x10ffff
            buf[i += 1] = UInt16(0xd7c0 + (c >>> 10))
            buf[i += 1] = UInt16(0xdc00 + (c & 0x3ff))
        else
            throw(ArgumentError("invalid Unicode character (0x$(hex(c)) > 0x10ffff)"))
        end
    end
    buf[i += 1] = 0 # NULL termination
    UTF16String(buf)
end
=#
function encode16(s::AbstractString)
    buf = UInt16[]
    for ch in s
        c = reinterpret(UInt32, ch)
        if c < 0x10000
            push!(buf, UInt16(c))
        elseif c <= 0x10ffff
            push!(buf, UInt16(0xd7c0 + (c>>10)))
            push!(buf, UInt16(0xdc00 + (c & 0x3ff)))
        else
            throw(ArgumentError("invalid Unicode character (0x$(hex(c)) > 0x10ffff)"))
        end
    end
    push!(buf, 0) # NULL termination
    UTF16String(buf)
end
#=
"""
@brief      reencodes an ASCII string using UTF-16 encoding

@param[in]  UTF16String
@param[in]  str::Vector{UInt8}
@param[in]  ASCIIString

@return     ::UTF16String
@throws     ArgumentError
"""
=#
function reencode_string(::Type{UTF16String}, str::Vector{UInt8}, ::Type{ASCIIString})
    len::Int = sizeof(str)
    buf = Vector{UInt16}(len+1)
    for i = 1:len
        @inbounds buf[i] = str[i]
    end
    @inbounds buf[len+1] = 0 # NULL termination
    UTF16String(buf)
end

const UTF_NO_LONG_NULL = 1      # don't accept 0xc0 0x80 for '\0'
const UTF_NO_SURROGATES = 2     # don't accept surrogate pairs in UTF-8/UTF-32
const UTF_ACCEPT_LONG = 4       # accept long encodings (other than long null in UTF-8)

const UTF_LONG = 1              # Long encodings are present
const UTF_LATIN1 = 2            # characters in range 0x80-0xFF present
const UTF_UNICODE2 = 4          # characters in range 0x100-0x7ff present
const UTF_UNICODE3 = 8          # characters in range 0x800-0xd7ff, 0xe000-0xffff
const UTF_UNICODE4 = 16         # non-BMP characters present
const UTF_SURROGATE = 32        # surrogate pairs present

#=
@doc """
@brief      Validates and calculates number of characters in UTF-8 string

@param[in]  UTF8String
@param[in]  str UTF-8 string
@param[in]  options flags to determine error handling (default 0)

@return     (total characters, 2-byte, 3-byte, 4-byte, flags)
@throws     ArgumentError
""" ->
=#
function check_string(::Type{UTF8String}, str::Vector{UInt8}, options::Integer=0)
    ch::UInt32 = 0
    c1::UInt32 = 0
    c2::UInt32 = 0
    c3::UInt32 = 0
    cntT::UInt = 0      # total # of characters
    cnt2::UInt = 0      # number of characters in the range (0x100-0x7ff)
    cnt3::UInt = 0      # number of characters in the range (0x800-0xd7ff,0xe000-0xffff)
    cnt4::UInt = 0      # number of characters in the range (0x10000-0x10ffff)
    flags::UInt = 0
    pos::UInt = 0
    len::UInt = sizeof(str)
    @inbounds while pos < len
        ch = str[pos += 1]
        cntT += 1
        if ch > 0x7f
            if ch < 0xe0
                (pos == len) && utf_errfunc(UTF_ERR_SHORT, pos, ch)
                c1 = str[pos += 1]
                !is_valid_continuation(c1) && utf_errfunc(UTF_ERR_CONT, pos, ch)
                ch = ((ch & 0x3f) << 6) | (c1 & 0x3f)
                if ch > 0x7f
                    cnt2 += 1
                    flags |= (ch > 0xff) ? UTF_UNICODE2 : UTF_LATIN1
                elseif (options & UTF_ACCEPT_LONG) != 0
                    flags |= UTF_LONG
                elseif (ch == 0) && ((options & UTF_NO_LONG_NULL) == 0)
                    flags |= UTF_LONG
                else
                    utf_errfunc(UTF_ERR_LONG, pos, ch)
                end
            elseif ch < 0xf0
                (pos + 2 > len) && utf_errfunc(UTF_ERR_SHORT, pos, ch)
                c1 = str[pos += 1]
                !is_valid_continuation(c1) && utf_errfunc(UTF_ERR_CONT, pos, c1)
                c2 = str[pos += 1]
                !is_valid_continuation(c2) && utf_errfunc(UTF_ERR_CONT, pos, c2)
                ch = ((ch & 0x0f) << 12) | ((c1 & 0x3f) << 6) | (c2 & 0x3f)
                # check for surrogate pairs, make sure correct
                if is_surrogate_char(ch)
                    !is_surrogate_lead(ch) && utf_errfunc(UTF_ERR_NOT_LEAD, pos-2, ch)
                    # next character *must* be a trailing surrogate character
                    (pos + 3 > len) && utf_errfunc(UTF_ERR_MISSING_SURROGATE, pos-2, ch)
                    c1 = str[pos += 1]
                    (c1 != 0xed) && utf_errfunc(UTF_ERR_NOT_TRAIL, pos, c1)
                    c2 = str[pos += 1]
                    !is_valid_continuation(c2) && utf_errfunc(UTF_ERR_CONT, pos, c2)
                    c3 = str[pos += 1]
                    !is_valid_continuation(c3) && utf_errfunc(UTF_ERR_CONT, pos, c3)
                    c3 = ((c1 & 0x0f) << 12) | ((c2 & 0x3f) << 6) | (c3 & 0x3f)
                    !is_surrogate_trail(c3) && utf_errfunc(UTF_ERR_NOT_TRAIL, pos-2, c3)
                    (options & UTF_NO_SURROGATES) != 0 && utf_errfunc(UTF_ERR_SURROGATE, pos-2, c3)
                    flags |= UTF_SURROGATE
                    cnt4 += 1
                elseif ch > 0x07ff
                    cnt3 += 1
                elseif (options & UTF_ACCEPT_LONG) != 0
                    flags |= UTF_LONG
                    cnt2 += 1
                else
                    utf_errfunc(UTF_ERR_LONG, pos-2, ch)
                end
            elseif ch < 0xf5
                (pos + 3 > len) && utf_errfunc(UTF_ERR_SHORT, pos, ch)
                c1 = str[pos += 1]
                !is_valid_continuation(c1) && utf_errfunc(UTF_ERR_CONT, pos, c1)
                c2 = str[pos += 1]
                !is_valid_continuation(c2) && utf_errfunc(UTF_ERR_CONT, pos, c2)
                c3 = str[pos += 1]
                !is_valid_continuation(c3) && utf_errfunc(UTF_ERR_CONT, pos, c3)
                ch = ((ch & 0x07) << 18) | ((c1 & 0x3f) << 12) | ((c2 & 0x3f) << 6) | (c3 & 0x3f)
                if ch > 0x10ffff
                    utf_errfunc(UTF_ERR_INVALID, pos-3, ch)
                elseif ch > 0xffff
                    cnt4 += 1
                elseif is_surrogate_char(ch)
                    utf_errfunc(UTF_ERR_SURROGATE, pos-3, ch)
                elseif (options & UTF_ACCEPT_LONG) != 0
                    flags |= UTF_LONG
                    if ch > 0x7ff
                        cnt3 += 1
                    elseif ch > 0x7f
                        cnt2 += 1
                    end
                else
                    utf_errfunc(UTF_ERR_LONG, pos-2, ch)
                end
            else
                utf_errfunc(UTF_ERR_INVALID, pos, ch)
            end
        end
    end
    cntT, cnt2, cnt3, cnt4, flags | (cnt3 == 0 ? 0 : UTF_UNICODE3) | (cnt4 == 0 ? 0 : UTF_UNICODE4)
end

#=
@doc """
@brief      Validates and calculates number of characters in UTF-16 string

@param[in]  UTF16String
@param[in]  str::UTF16String
@param[in]  options flags to determine error handling (default 0)

@return     (total characters, 2-byte, 3-byte, 4-byte, flags)
@throws     ArgumentError
""" ->
=#
function check_string(::Type{UTF16String}, str::Vector{UInt16}, options::Integer=0)
    ch::UInt32 = 0
    cntT::UInt = 0      # total # of characters
    cnt2::UInt = 0      # number of characters in the range (0x80-0x7ff)
    cnt3::UInt = 0      # number of characters in the range (0x800-0xd7ff,0xe000-0xffff)
    cnt4::UInt = 0      # number of characters in the range (0x10000-0x10ffff)
    flags::UInt = 0
    pos::UInt = 0
    len::UInt = sizeof(str) >>> 1
    @inbounds begin
    str[len] == 0 && (len -= 1)
    while pos < len
        ch = str[pos += 1]
        cntT += 1
        if ch > 0x7f
            if ch < 0x100
                cnt2 += 1
                flags |= UTF_LATIN1
            elseif ch < 0x800
                cnt2 += 1
                flags |= UTF_UNICODE2
            elseif !is_surrogate_char(ch)
                cnt3 += 1
            elseif is_surrogate_lead(ch)
                pos == len && utf_errfunc(UTF_ERR_MISSING_SURROGATE, pos, ch)
                # next character *must* be a trailing surrogate character
                ch = str[pos += 1]
                !is_surrogate_trail(ch) && utf_errfunc(UTF_ERR_NOT_TRAIL, pos, ch)
                cnt4 += 1
            else
                utf_errfunc(UTF_ERR_NOT_LEAD, pos, ch)
            end
        end
    end
    end
    cntT, cnt2, cnt3, cnt4, flags | (cnt3 == 0 ? 0 : UTF_UNICODE3) | (cnt4 == 0 ? 0 : UTF_UNICODE4)
end

#=
"""
@brief      Validates and calculates number of characters in UTF-32 string

@param[in]  UTF32String
@param[in]  str::Vector{UInt32}
@param[in]  options flags to determine error handling (default 0)

@return     (total characters, 2-byte, 3-byte, 4-byte, flags)
@throws     ArgumentError
"""
=#
function check_string(::Type{UTF32String}, str::Vector{UInt32}, options::Integer=0)
    cntT::UInt = 0      # total # of characters
    cnt2::UInt = 0      # number of characters in the range (0x80-0x7ff)
    cnt3::UInt = 0      # number of characters in the range (0x800-0xd7ff,0xe000-0xffff)
    cnt4::UInt = 0      # number of characters in the range (0x10000-0x10ffff)
    pos::UInt = 0
    flags::UInt = 0
    len::UInt = (sizeof(str) >>> 2)
    @inbounds begin
    str[len] == 0 && (len -= 1)
    while pos < len
        ch = str[pos += 1]
        cntT += 1
        if ch > 0x7f
            if ch < 0x100
                cnt2 += 1
                flags |= UTF_LATIN1
            elseif ch < 0x800
                cnt2 += 1
                flags |= UTF_UNICODE2
            elseif ch > 0xffff
                (ch > 0x10ffff) && utf_errfunc(UTF_ERR_INVALID, pos, ch)
                cnt4 += 1
            elseif !is_surrogate_char(ch)
                cnt3 += 1
            elseif is_surrogate_lead(ch)
                pos == len && utf_errfunc(UTF_ERR_MISSING_SURROGATE, pos, ch)
                # next character *must* be a trailing surrogate character
                ch = str[pos += 1]
                !is_surrogate_trail(ch) && utf_errfunc(UTF_ERR_NOT_TRAIL, pos, ch)
                cnt4 += 1
                (options & UTF_NO_SURROGATES) != 0 && utf_errfunc(UTF_ERR_SURROGATE, pos, ch)
                flags |= UTF_SURROGATE
            else
                utf_errfunc(UTF_ERR_NOT_LEAD, pos, ch)
            end
        end
    end
    end
    cntT, cnt2, cnt3, cnt4, flags | (cnt3 == 0 ? 0 : UTF_UNICODE3) | (cnt4 == 0 ? 0 : UTF_UNICODE4)
end

#=
"""
@brief      Validates and calculates number of characters in UTF-32 string

@param[in]  UTF32String
@param[in]  str::AbstractString
@param[in]  options flags to determine error handling (default 0)

@return     (total characters, 2-byte, 3-byte, 4-byte, flags)
@throws     ArgumentError
"""
=#
function check_string(::Type{UTF32String}, str::AbstractString, options::Integer=0)
    cntT::UInt = 0      # total # of characters
    cnt2::UInt = 0      # number of characters in the range (0x80-0x7ff)
    cnt3::UInt = 0      # number of characters in the range (0x800-0xd7ff,0xe000-0xffff)
    cnt4::UInt = 0      # number of characters in the range (0x10000-0x10ffff)
    pos::UInt = start(str)
    len::UInt = endof(str)
    flags::UInt = 0
    while pos < len
        cntT += 1
        ch, pos = next(str, pos)
        if ch > 0x7f
            if ch < 0x100
                cnt2 += 1
                flags |= UTF_LATIN1
            elseif ch < 0x800
                cnt2 += 1
                flags |= UTF_UNICODE2
            elseif ch > 0xffff
                (ch > 0x10ffff) && utf_errfunc(UTF_ERR_INVALID, pos, ch)
                cnt4 += 1
            elseif !is_surrogate_char(ch)
                cnt3 += 1
            elseif is_surrogate_lead(ch)
                pos == len && utf_errfunc(UTF_ERR_MISSING_SURROGATE, pos, ch)
                # next character *must* be a trailing surrogate character
                ch, pos = next(str, pos)
                !is_surrogate_trail(ch) && utf_errfunc(UTF_ERR_NOT_TRAIL, pos, ch)
                cnt4 += 1
                (options & UTF_NO_SURROGATES) != 0 && utf_errfunc(UTF_ERR_SURROGATE, pos, ch)
                flags |= UTF_SURROGATE
            else
                utf_errfunc(UTF_ERR_NOT_LEAD, pos, ch)
            end
        end
    end
    cntT, cnt2, cnt3, cnt4, flags | (cnt3 == 0 ? 0 : UTF_UNICODE3) | (cnt4 == 0 ? 0 : UTF_UNICODE4)
end

# const errFunc = cfunction(Base.utf_errfunc, Cssize_t, (Int, Csize_t, Cuint))

#=
@doc """
@brief      Converts a UTF-8 encoded string to UTF-16 encoding

@param[in]  UTF16String
@param[in]  str::Vector{UInt8}
@param[in]  UTF8String

@return     ::UTF16String
@throws     ArgumentError
""" ->
=#
function reencode_string(::Type{UTF16String}, str::Vector{UInt8}, ::Type{UTF8String})
    len::UInt = sizeof(str)
    # handle zero length string quickly
    if len == 0
        buf = Vector{UInt16}(1)
        buf[1] = 0
        return UTF16String(buf)
    end
    # Check that is correct UTF-8 encoding and get number of words needed
    len, cnt2::UInt, cnt3::UInt, cnt4::UInt, flags::UInt = check_string(UTF8String, str)
    # Optimize case where no characters > 0x7f
    flags == 0 && return reencode_string(UTF16String, str, ASCIIString)
    len += cnt4
    buf = Vector{UInt16}(len+1)
    @inbounds begin
    buf[len+1] = 0
    out::UInt = 0
    pos::UInt = 0
    while out < len
        ch::UInt32 = str[pos += 1]
        if ch > 0x7f
            if ch < 0xe0
                ch = ((ch & 0x1f) << 6) | (str[pos += 1] & 0x3f)
            elseif ch < 0xf0
                ch = ((ch & 0xf) << 12) | (UInt32(str[pos+1] & 0x3f) << 6) | (str[pos+2] & 0x3f)
                pos += 2
            else
                ch = ((ch & 0x7) << 18) | (UInt32(str[pos+1] & 0x3f) << 12) | (UInt32(str[pos+2] & 0x3f) << 6) | (str[pos+3] & 0x3f)
                pos += 3
                buf[out += 1] = 0xd7c0 + (ch >>> 10)
                ch = 0xdc00 + (ch & 0x3ff)
            end
        end
        buf[out += 1] = ch
    end
    end
    UTF16String(buf)
end

#=
@doc """
@brief      Reencodes a UTF-16 encoded string using UTF-8 encoding

@param[in]  UTF8String
@param[in]  str::Vector{UInt16}
@param[in]  UTF16String

@return     ::UTF8String
@throws     ArgumentError
""" ->
=#
function reencode_string(::Type{UTF8String}, str::Vector{UInt16}, ::Type{UTF16String})
    # handle zero length string quickly
    sizeof(str) <= 2 && return UTF8String("")
    # get number of bytes to allocate
    len::UInt, cnt2::UInt, cnt3::UInt, cnt4::UInt, flags::UInt = check_string(UTF16String, str)
    len += cnt2 + cnt3*2 + cnt4*3
    buf = Vector{UInt8}(len)
    @inbounds begin
    if flags == 0
        # If copy! is fixed to handle narrowing operations efficiently, change this to use copy!
        # 7-bit ASCII characters only! Do optimized conversion
        for pos=1:len
            buf[pos] = str[pos]
        end
    else
        out::UInt = 0
        pos::UInt = 0
        while out < len
            ch::UInt32 = str[pos += 1]
            if ch > 0x7f
                if ch < 0x800
                    buf[out += 1] = 0xc0 | (ch >>> 6)
                elseif is_surrogate_char(ch)
                    ch = 0x10000 + (((ch & 0x3ff) << 10) | (str[pos += 1] & 0x3ff))
                    buf[out += 1] = 0xf0 | (ch >>> 18)
                    buf[out += 1] = 0x80 | ((ch >>> 12) & 0x3f)
                    buf[out += 1] = 0x80 | ((ch >>> 6) & 0x3f)
                else
                    buf[out += 1] = 0xe0 | ((ch >>> 12) & 0x3f)
                    buf[out += 1] = 0x80 | ((ch >>> 6) & 0x3f)
                end
                ch = 0x80 | (ch & 0x3f)
            end
            buf[out += 1] = ch
        end
    end
    end
    UTF8String(buf)
end

utf16(x) = convert(UTF16String, x)
convert(::Type{UTF16String}, s::UTF16String) = s
convert(::Type{UTF16String}, s::AbstractString) = encode16(s)
convert(::Type{Vector{UInt16}}, s::UTF16String) = s.data
convert(::Type{Array{UInt16}}, s::UTF16String) = s.data

convert(::Type{UTF8String}, str::Vector{UInt16}) = reencode_string(UTF8String, str, UTF16String)
convert(::Type{UTF8String},  str::UTF16String) = reencode_string(UTF8String, str.data, UTF16String)
convert(::Type{UTF16String}, str::ASCIIString) = reencode_string(UTF16String, str.data, ASCIIString)
convert(::Type{UTF16String}, str::UTF8String) = reencode_string(UTF16String, str.data, UTF8String)

sizeof(s::UTF16String) = sizeof(s.data) - sizeof(UInt16)
unsafe_convert{T<:Union(Int16,UInt16)}(::Type{Ptr{T}}, s::UTF16String) =
    convert(Ptr{T}, pointer(s))

function is_valid_utf16(data::AbstractArray{UInt16})
    i = 1
    n = length(data) # this may include NULL termination; that's okay
    while i < n # check for unpaired surrogates
        if utf16_is_lead(data[i]) && utf16_is_trail(data[i+1])
            i += 2
        elseif utf16_is_surrogate(data[i])
            return false
        else
            i += 1
        end
    end
    return i > n || !utf16_is_surrogate(data[i])
end

is_valid_utf16(s::UTF16String) = is_valid_utf16(s.data)

function convert(::Type{UTF16String}, data::AbstractVector{UInt16})
    !is_valid_utf16(data) && throw(ArgumentError("invalid UTF16 data"))
    len = length(data)
    d = Array(UInt16, len + 1)
    d[end] = 0 # NULL terminate
    UTF16String(copy!(d,1, data,1, len))
end

convert(T::Type{UTF16String}, data::AbstractArray{UInt16}) =
    convert(T, reshape(data, length(data)))

convert(T::Type{UTF16String}, data::AbstractArray{Int16}) =
    convert(T, reinterpret(UInt16, data))

function convert(T::Type{UTF16String}, bytes::AbstractArray{UInt8})
    isempty(bytes) && return UTF16String(UInt16[0])
    isodd(length(bytes)) && throw(ArgumentError("odd number of bytes"))
    data = reinterpret(UInt16, bytes)
    # check for byte-order mark (BOM):
    if data[1] == 0xfeff        # native byte order
        d = Array(UInt16, length(data))
        copy!(d,1, data,2, length(data)-1)
    elseif data[1] == 0xfffe    # byte-swapped
        d = Array(UInt16, length(data))
        for i = 2:length(data)
            d[i-1] = bswap(data[i])
        end
    else
        d = Array(UInt16, length(data) + 1)
        copy!(d,1, data,1, length(data)) # assume native byte order
    end
    d[end] = 0 # NULL terminate
    !is_valid_utf16(d) && throw(ArgumentError("invalid UTF16 data"))
    UTF16String(d)
end

utf16(p::Ptr{UInt16}, len::Integer) = utf16(pointer_to_array(p, len))
utf16(p::Ptr{Int16}, len::Integer) = utf16(convert(Ptr{UInt16}, p), len)
function utf16(p::Union(Ptr{UInt16}, Ptr{Int16}))
    len = 0
    while unsafe_load(p, len+1) != 0; len += 1; end
    utf16(p, len)
end
