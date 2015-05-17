# This file is a part of Julia. License is MIT: http://julialang.org/license

next(s::UTF32String, i::Int) = (s.data[i], i+1)
endof(s::UTF32String) = length(s.data) - 1
length(s::UTF32String) = length(s.data) - 1

function utf32(c::Integer...)
    a = Array(Char, length(c) + 1)
    for i = 1:length(c)
        @inbounds a[i] = Char(c[i])
    end
    @inbounds a[end] = Char(0)
    UTF32String(a)
end

#=
"""
@brief      Reencodes an ASCII string using UTF-32 encoding

@param[in]  str::Vector{UInt8}

@return     ::UTF32String
@throws     ArgumentError
"""
=#
function reencode_string_a_32(str::Vector{UInt8})
    len::Int = sizeof(str)
    buf = Vector{Char}(len+1)
    for i = 1:len
        @inbounds buf[i] = Char(str[i])
    end
    @inbounds buf[len+1] = Char(0) # NULL termination
    UTF32String(buf)
end

#=
"""
@brief      Reencodes a UTF-8 encoded string using UTF-32 encoding

@param[in]  str::Vector{UInt8}

@return     ::UTF32String
@throws     ArgumentError
"""
=#
function reencode_string_8_32(str::Vector{UInt8})
    len::UInt = sizeof(str)
    # handle zero length string quickly
    if len == 0
        buf = Vector{Char}(1)
        @inbounds buf[1] = Char(0)
        return UTF32String(buf)
    end
    # Validate UTF-8 encoding, and get number of words to create
    len, cnt2::UInt, cnt3::UInt, cnt4::UInt, flags::UInt = check_string_utf8(str)
    buf = Vector{Char}(len+1)
    @inbounds begin
    buf[len+1] = Char(0) # NULL termination
    if flags == 0
        # 7-bit ASCII characters only! Do optimized conversion
        for pos=1:len
            buf[pos] = Char(str[pos])
        end
    else
        # has multi-byte UTF-8 sequences
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
		    if is_surrogate_lead(ch)
			ch = (((ch & 0x3ff) << 10)
				| (((UInt32(str[pos+1] & 0xf) << 12)
				|  (UInt32(str[pos+2] & 0x3f) << 6)
				|  (str[pos+3] & 0x3f)) & 0x3ff)) + 0x10000
			pos += 3
		    end
                else
                    ch = ((ch & 0x7) << 18) | (UInt32(str[pos+1] & 0x3f) << 12) | (UInt32(str[pos+2] & 0x3f) << 6) | (str[pos+3] & 0x3f)
                    pos += 3
                end
            end
            buf[out += 1] = Char(ch)
        end
    end
    end
    UTF32String(buf)
end

#=
"""
@brief      Reencodes a UTF-16 encoded string using UTF-32 encoding

@param[in]  str::Vector{UInt16}

@return     ::UTF32String
@throws     ArgumentError
"""
=#
function reencode_string_16_32(str::Vector{UInt16})
    len::UInt = (sizeof(str) >>> 1) - 1 # account for trailing \0
    # handle zero length string quickly
    if len == 0
        buf = Vector{Char}(1)
        @inbounds buf[1] = Char(0)
        return UTF32String(buf)
    end
    # get number of words to create
    len, cnt2::UInt, cnt3::UInt, cnt4::UInt, flags::UInt = check_string_utf16(str)
    buf = Vector{Char}(len+1)
    @inbounds begin
    buf[len+1] = Char(0) # NULL termination
    if (flags & UTF_UNICODE4) == 0
        # No surrogate pairs, do optimized copy
        for pos=1:len
            buf[pos] = Char(str[pos])
        end
    else
        # has multi-byte UTF-8 sequences
        out::UInt = 0
        pos::UInt = 0
        while out < len
            ch::UInt32 = str[pos += 1]
            # check for surrogate pair
            if is_surrogate_lead(ch)
                ch = 0x10000 + (((ch & 0x3ff) << 10) | (str[pos += 1] & 0x3ff))
            end
            buf[out += 1] = Char(ch)
        end
    end
    end
    UTF32String(buf)
end

#=
"""
@brief      Reencode a UTF-32 encoded string using UTF-16 encoding

@param[in]  str::Vector{Char}

@return     ::UTF16String
@throws     ArgumentError
"""
=#
function reencode_string_32_16(str::Vector{UInt32})
    len::UInt = (sizeof(str) >>> 2) - 1 # account for trailing \0
    # handle zero length string quickly
    if len == 0 ; buf = Vector{UInt16}(1) ; buf[1] = 0 ; return UTF16String(buf) ; end
    # get number of words to allocate
    len, cnt2::UInt, cnt3::UInt, cnt4::UInt, flags::UInt = check_string_utf32(str)
    buf = Vector{UInt16}(len+cnt4+1)
    @inbounds begin
    buf[len+cnt4+1] = 0 # NULL termination
    if cnt4 == 0        # optimized path, no surrogates
        for pos = 1:len
            buf[pos] = UInt16(str[pos])
        end
    else
        out::UInt = 0
        pos::UInt = 0
        while pos < len
            ch = UInt32(str[pos += 1])
            if ch > 0xffff
                buf[out += 1] = 0xd7c0 + (ch >>> 10)
                ch = 0xdc00 + (ch & 0x3ff)
            end
            buf[out += 1] = ch
        end
    end
    end
    UTF16String(buf)
end

#=
"""
@brief      reencode_string a UTF-32 encoded string using UTF-8 encoding

@param[in]  str::Union(UTF32String, Vector{Char})

@return     ::UTF8String
@throws     ArgumentError
"""
=#
function reencode_string_32_8(str::Vector{UInt32})
    len::UInt = (sizeof(str) >>> 2) - 1 # account for trailing \0
    # handle zero length string quickly
    len == 0 && return UTF8String("")
    # get number of bytes to allocate
    len, cnt2::UInt, cnt3::UInt, cnt4::UInt, flags::UInt = check_string_utf32(str)
    len += cnt2 + cnt3*2 + cnt4*3
    buf = Vector{UInt8}(len)
    @inbounds begin
    if flags == 0
        # If copy! is fixed to handle narrowing operations efficiently, change this to use copy!
        # 7-bit ASCII characters only! Do optimized conversion
        for pos=1:len
             buf[pos] = reinterpret(UInt32, str[pos])
        end
    else
        out::UInt = 0
        pos::UInt = 0
        while out < len
            ch::UInt32 = reinterpret(UInt32, str[pos += 1])
            if ch > 0x7f
                if ch < 0x800
                    buf[out += 1] = 0xc0 | (ch >>> 6)
                elseif is_surrogate_char(ch)
                    ch = 0x10000 + (((ch & 0x3ff) << 10) | (reinterpret(UInt32,str[pos += 1]) & 0x3ff))
                    buf[out += 1] = 0xf0 | (ch >>> 18)
                    buf[out += 1] = 0x80 | ((ch >>> 12) & 0x3f)
                    buf[out += 1] = 0x80 | ((ch >>> 6) & 0x3f)
                elseif ch > 0xffff
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

utf32(x) = convert(UTF32String, x)

convert(::Type{UTF32String}, c::Char)             = UTF32String(Char[c, Char(0)])
convert(::Type{UTF32String}, str::UTF32String)    = str

convert(::Type{UTF8String},  str::UTF32String)    = reencode_string_32_8(reinterpret(UInt32, str.data))
convert(::Type{UTF8String},  str::Vector{Char})   = reencode_string_32_8(reinterpret(UInt32, str))
convert(::Type{UTF8String},  str::Vector{UInt32}) = reencode_string_32_8(str)

convert(::Type{UTF16String}, str::UTF32String)    = reencode_string_32_16(reinterpret(UInt32, str.data))
convert(::Type{UTF16String}, str::Vector{Char})   = reencode_string_32_16(reinterpret(UInt32, str))
convert(::Type{UTF16String}, str::Vector{UInt32}) = reencode_string_32_16(str)

convert(::Type{UTF32String}, str::ASCIIString)    = reencode_string_a_32(str.data)
convert(::Type{UTF32String}, str::UTF8String)     = reencode_string_8_32(str.data)
convert(::Type{UTF32String}, str::UTF16String)    = reencode_string_16_32(str.data)


# TODO, this doesn't check for surrogate pairs, which can be present
# in poorly encoded data, and doesn't check validity of Chars
function convert(::Type{UTF32String}, s::AbstractString)
    a = Array(Char, length(s) + 1)
    i = 0
    for c in s
        a[i += 1] = c
    end
    a[end] = Char(0) # NULL terminate
    UTF32String(a)
end

# TODO, this doesn't check for surrogate pairs, which can be present
# in poorly encoded data, and doesn't check validity of Chars
function convert(::Type{UTF32String}, data::AbstractVector{Char})
    len = length(data)
    d = Array(Char, len + 1)
    d[end] = Char(0) # NULL terminate
    UTF32String(copy!(d,1, data,1, len))
end

convert{T<:Union(Int32,UInt32)}(::Type{UTF32String}, data::AbstractVector{T}) =
    convert(UTF32String, reinterpret(Char, data))

convert{T<:AbstractString}(::Type{T}, v::AbstractVector{Char}) = convert(T, utf32(v))

# specialize for performance reasons:
function convert{T<:ByteString}(::Type{T}, data::AbstractVector{Char})
    s = IOBuffer(Array(UInt8,length(data)), true, true)
    truncate(s,0)
    for x in data
        print(s, x)
    end
    convert(T, takebuf_string(s))
end

convert(::Type{Array{Char,1}}, s::UTF32String) = s.data
convert(::Type{Array{Char}}, s::UTF32String) = s.data

reverse(s::UTF32String) = UTF32String(reverse!(copy(s.data), 1, length(s)))

sizeof(s::UTF32String) = sizeof(s.data) - sizeof(Char)
unsafe_convert{T<:Union(Int32,UInt32,Char)}(::Type{Ptr{T}}, s::UTF32String) =
    convert(Ptr{T}, pointer(s))

# TODO, this doesn't check for surrogate pairs, which can be present
# in poorly encoded data, and doesn't check validity of input data
function convert(T::Type{UTF32String}, bytes::AbstractArray{UInt8})
    isempty(bytes) && return UTF32String(Char[0])
    length(bytes) & 3 != 0 && throw(ArgumentError("need multiple of 4 bytes"))
    data = reinterpret(Char, bytes)
    # check for byte-order mark (BOM):
    if data[1] == Char(0x0000feff) # native byte order
        d = Array(Char, length(data))
        copy!(d,1, data, 2, length(data)-1)
    elseif data[1] == Char(0xfffe0000) # byte-swapped
        d = Array(Char, length(data))
        for i = 2:length(data)
            d[i-1] = bswap(data[i])
        end
    else
        d = Array(Char, length(data) + 1)
        copy!(d, 1, data, 1, length(data)) # assume native byte order
    end
    d[end] = Char(0) # NULL terminate
    UTF32String(d)
end

function isvalid(::Type{UTF32String}, str::Union(Vector{Char}, Vector{UInt32}))
    for i=1:length(str)
        @inbounds if !isvalid(Char, reinterpret(UInt32, str[i])) ; return false ; end
    end
    return true
end
isvalid(str::Vector{Char}) = isvalid(UTF32String, str)
isvalid{T<:Union(ASCIIString,UTF8String,UTF16String,UTF32String)}(str::T) = isvalid(T, str.data)

utf32(p::Ptr{Char}, len::Integer) = utf32(pointer_to_array(p, len))
utf32(p::Union(Ptr{UInt32}, Ptr{Int32}), len::Integer) = utf32(convert(Ptr{Char}, p), len)
function utf32(p::Union(Ptr{Char}, Ptr{UInt32}, Ptr{Int32}))
    len = 0
    while unsafe_load(p, len+1) != 0; len += 1; end
    utf32(p, len)
end

function map(f, s::UTF32String)
    d = s.data
    out = similar(d)
    out[end] = Char(0)

    for i = 1:(length(d)-1)
        c2 = f(d[i])
        if !isa(c2, Char)
            throw(ArgumentError("map(f,s::AbstractString) requires f to return Char; try map(f,collect(s)) or a comprehension instead"))
        end
        out[i] = (c2::Char)
    end
    UTF32String(out)
end
