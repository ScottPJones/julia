# This file is a part of Julia. License is MIT: https://julialang.org/license

"""
    AbstractDateToken

A token used in parsing or formatting a date time string. Each subtype must
define the tryparsenext and format methods.

"""
abstract type AbstractDateToken end

"""
    tryparsenext(tok::AbstractDateToken, str::String, i::Int, len::Int, locale::DateLocale)

`tryparsenext` parses for the `tok` token in `str` starting at index `i`.
`len` is the length of the string.  parsing can be optionally based on the
`locale`. If a `tryparsenext` method does not need a locale, it can leave
the argument out in the method definition.

If parsing succeeds, returns a tuple of 2 elements `(res, idx)`, where:

* `res` is the result of the parsing.
* `idx::Int`, is the index _after_ the index at which parsing ended.
"""
function tryparsenext end

"""
    format(io::IO, tok::AbstractDateToken, dt::TimeType, locale)

Format the `tok` token from `dt` and write it to `io`. The formatting can
be based on `locale`.

All subtypes of `AbstractDateToken` must define this method in order
to be able to print a Date / DateTime object according to a `DateFormat`
containing that token.
"""
format(io::IO, tok::AbstractDateToken, dt::TimeType, locale)

# fallback to tryparsenext/format methods that don't care about locale
@inline function tryparsenext(d::AbstractDateToken, str, i, len, locale)
    return tryparsenext(d, str, i, len)
end

function Base.string(t::Time)
    h, mi, s = hour(t), minute(t), second(t)
    hh = lpad(h, 2, "0")
    mii = lpad(mi, 2, "0")
    ss = lpad(s, 2, "0")
    nss = tons(Millisecond(t)) + tons(Microsecond(t)) + tons(Nanosecond(t))
    ns = nss == 0 ? "" : rstrip(@sprintf("%.9f", nss / 1e+9)[2:end], '0')
    return "$hh:$mii:$ss$ns"
end

Base.show(io::IO, ::MIME"text/plain", t::Time) = print(io, t)
Base.print(io::IO, t::Time) = print(io, string(t))

function Base.show(io::IO, t::Time)
    if get(io, :compact, false)::Bool
        print(io, t)
    else
        values = [
            hour(t)
            minute(t)
            second(t)
            millisecond(t)
            microsecond(t)
            nanosecond(t)
        ]
        index = something(findlast(!iszero, values), 1)

        print(io, Time, "(")
        for i in 1:index
            show(io, values[i])
            i != index && print(io, ", ")
        end
        print(io, ")")
    end
end

@inline function format(io, d::AbstractDateToken, dt, locale)
    format(io, d, dt)
end

# Information for parsing and formatting date time values.
struct DateFormat{S, T<:Tuple}
    tokens::T
    locale::DateLocale
end

### Token types ###

struct DatePart{letter} <: AbstractDateToken
    width::Int
    fixed::Bool
end

@inline min_width(d::DatePart) = d.fixed ? d.width : 1
@inline max_width(d::DatePart) = d.fixed ? d.width : 0

function _show_content(io::IO, d::DatePart{c}) where c
    for i = 1:d.width
        print(io, c)
    end
end

function Base.show(io::IO, d::DatePart{c}) where c
    print(io, "DatePart(")
    _show_content(io, d)
    print(io, ")")
end

### Parse tokens

for c in "yY"
    @eval begin
        @inline function tryparsenext(d::DatePart{$c}, str, i, len)
            val = tryparsenext_sign(str, i, len)
            if val !== nothing
                coefficient, i = val
            else
                coefficient = 1
            end
            # The sign character does not affect fixed length `DatePart`s
            val = tryparsenext_base10(str, i, len, min_width(d), max_width(d))
            val === nothing && return nothing
            y, ii = val
            return y * coefficient, ii
        end
    end
end

for c in "mdHIMS"
    @eval begin
        @inline function tryparsenext(d::DatePart{$c}, str, i, len)
            return tryparsenext_base10(str, i, len, min_width(d), max_width(d))
        end
    end
end

function tryparsenext(d::DatePart{'p'}, str, i, len)
    i+1 > len && return nothing
    c, ii = iterate(str, i)::Tuple{Char, Int}
    ap = lowercase(c)
    (ap == 'a' || ap == 'p') || return nothing
    c, ii = iterate(str, ii)::Tuple{Char, Int}
    lowercase(c) == 'm' || return nothing
    return ap == 'a' ? AM : PM, ii
end

for (tok, fn) in zip("uUeE", Any[monthabbr_to_value, monthname_to_value, dayabbr_to_value, dayname_to_value])
    @eval @inline function tryparsenext(d::DatePart{$tok}, str, i, len, locale)
        next = tryparsenext_word(str, i, len, locale, max_width(d))
        next === nothing && return nothing
        word, i = next
        val = $fn(word, locale)
        val == 0 && return nothing
        return val, i
    end
end

# 3-digit (base 10) number following a decimal point. For InexactError below.
struct Decimal3 end

@inline function tryparsenext(d::DatePart{'s'}, str, i, len)
    val = tryparsenext_base10(str, i, len, min_width(d), max_width(d))
    val === nothing && return nothing
    ms0, ii = val
    len = ii - i
    if len > 3
        ms, r = divrem(ms0, Int64(10) ^ (len - 3))
        r == 0 || return nothing
    else
        ms = ms0 * Int64(10) ^ (3 - len)
    end
    return ms, ii
end

### Format tokens

hour12(dt) = let h = hour(dt); h > 12 ? h - 12 : h == 0 ? 12 : h; end

for (c, fn) in zip("YmdHIMS", Any[year, month, day, hour, hour12, minute, second])
    @eval function format(io, d::DatePart{$c}, dt)
        print(io, string($fn(dt), base = 10, pad = d.width))
    end
end

for (tok, fn) in zip("uU", Any[monthabbr, monthname])
    @eval function format(io, d::DatePart{$tok}, dt, locale)
        print(io, $fn(month(dt), locale))
    end
end

function format(io, d::DatePart{'p'}, dt, locale)
    ampm = hour(dt) < 12 ? "AM" : "PM" # fixme: locale-specific?
    print(io, ampm)
end

for (tok, fn) in zip("eE", Any[dayabbr, dayname])
    @eval function format(io, ::DatePart{$tok}, dt, locale)
        print(io, $fn(dayofweek(dt), locale))
    end
end

@inline function format(io, d::DatePart{'y'}, dt)
    y = year(dt)
    n = d.width

    # the last n digits of y
    # will be 0 padded if y has less than n digits
    str = string(y, base = 10, pad = n)
    l = lastindex(str)
    if l == n
        # fast path
        print(io, str)
    else
        print(io, SubString(str, l - (n - 1), l))
    end
end

function format(io, d::DatePart{'s'}, dt)
    ms = millisecond(dt)
    if ms % 100 == 0
        str = string(div(ms, 100))
    elseif ms % 10 == 0
        str = string(div(ms, 10), pad = 2)
    else
        str = string(ms, pad = 3)
    end

    print(io, rpad(str, d.width, '0'))
end

### Delimiters

struct Delim{T, length} <: AbstractDateToken
    d::T
end

Delim(d::T) where {T<:AbstractChar} = Delim{T, 1}(d)
Delim(d::String) = Delim{String, length(d)}(d)

@inline function tryparsenext(d::Delim{<:AbstractChar, N}, str, i::Int, len) where N
    for j = 1:N
        i > len && return nothing
        next = iterate(str, i)
        @assert next !== nothing
        c, i = next
        c != d.d && return nothing
    end
    return true, i
end

@inline function tryparsenext(d::Delim{String, N}, str, i::Int, len) where N
    i1 = i
    i2 = firstindex(d.d)
    for j = 1:N
        if i1 > len
            return nothing
        end
        next1 = iterate(str, i1)
        @assert next1 !== nothing
        c1, i1 = next1
        next2 = iterate(d.d, i2)
        @assert next2 !== nothing
        c2, i2 = next2
        if c1 != c2
            return nothing
        end
    end
    return true, i1
end

@inline function format(io, d::Delim, dt, locale)
    print(io, d.d)
end

function _show_content(io::IO, d::Delim{<:AbstractChar, N}) where N
    if d.d in keys(CONVERSION_SPECIFIERS)
        for i = 1:N
            print(io, '\\', d.d)
        end
    else
        for i = 1:N
            print(io, d.d)
        end
    end
end

function _show_content(io::IO, d::Delim)
    for c in d.d
        if c in keys(CONVERSION_SPECIFIERS)
            print(io, '\\')
        end
        print(io, c)
    end
end

function Base.show(io::IO, d::Delim)
    print(io, "Delim(")
    _show_content(io, d)
    print(io, ")")
end

### DateFormat construction

abstract type DayOfWeekToken end # special addition to Period types

# Map conversion specifiers or character codes to tokens.
# Note: Allow addition of new character codes added by packages
const CONVERSION_SPECIFIERS = Dict{Char, Type}(
    'y' => Year,
    'Y' => Year,
    'm' => Month,
    'u' => Month,
    'U' => Month,
    'e' => DayOfWeekToken,
    'E' => DayOfWeekToken,
    'd' => Day,
    'H' => Hour,
    'I' => Hour,
    'M' => Minute,
    'S' => Second,
    's' => Millisecond,
    'p' => AMPM,
)

# Default values are needed when a conversion specifier is used in a DateFormat for parsing
# and we have reached the end of the input string.
# Note: Allow `Any` value as a default to support extensibility
const CONVERSION_DEFAULTS = IdDict{Type, Any}(
    Year => Int64(1),
    Month => Int64(1),
    DayOfWeekToken => Int64(0),
    Day => Int64(1),
    Hour => Int64(0),
    Minute => Int64(0),
    Second => Int64(0),
    Millisecond => Int64(0),
    Microsecond => Int64(0),
    Nanosecond => Int64(0),
    AMPM => TWENTYFOURHOUR,
)

# Specifies the required fields in order to parse a TimeType
# Note: Allows for addition of new TimeTypes
const CONVERSION_TRANSLATIONS = IdDict{Type, Any}(
    Date => (Year, Month, Day),
    DateTime => (Year, Month, Day, Hour, Minute, Second, Millisecond, AMPM),
    Time => (Hour, Minute, Second, Millisecond, Microsecond, Nanosecond, AMPM),
)

# The `DateFormat(format, locale)` method just below consumes the following Regex.
# Constructing this Regex is fairly expensive; doing so in the method itself can
# consume half or better of `DateFormat(format, locale)`'s runtime. So instead we
# construct and cache it outside the method body. Note, however, that when
# `keys(CONVERSION_SPECIFIERS)` changes, the cached Regex must be updated accordingly;
# hence the mutability (Ref-ness) of the cache, the helper method with which to populate
# the cache, the cache of the hash of `keys(CONVERSION_SPECIFIERS)` (to facilitate checking
# for changes), and the lock (to maintain consistency of these objects across threads when
# threads simultaneously modify `CONVERSION_SPECIFIERS` and construct `DateFormat`s).
function compute_dateformat_regex(conversion_specifiers)
    letters = String(collect(keys(conversion_specifiers)))
    return Regex("(?<!\\\\)([\\Q$letters\\E])\\1*")
end
const DATEFORMAT_REGEX_LOCK = ReentrantLock()
const DATEFORMAT_REGEX_HASH = Ref(hash(keys(CONVERSION_SPECIFIERS)))
const DATEFORMAT_REGEX_CACHE = Ref(compute_dateformat_regex(CONVERSION_SPECIFIERS))

"""
    DateFormat(format::AbstractString, locale="english")

Construct a date formatting object that can be used for parsing date strings or
formatting a date object as a string. The following character codes can be used to construct the `format`
string:

| Code       | Matches   | Comment                                                       |
|:-----------|:----------|:--------------------------------------------------------------|
| `Y`        | 1996, 96  | Returns year of 1996, 0096                                    |
| `y`        | 1996, 96  | Same as `Y` on `parse` but discards excess digits on `format` |
| `m`        | 1, 01     | Matches 1 or 2-digit months                                   |
| `u`        | Jan       | Matches abbreviated months according to the `locale` keyword  |
| `U`        | January   | Matches full month names according to the `locale` keyword    |
| `d`        | 1, 01     | Matches 1 or 2-digit days                                     |
| `H`        | 00        | Matches hours (24-hour clock)                                 |
| `I`        | 00        | For outputting hours with 12-hour clock                       |
| `M`        | 00        | Matches minutes                                               |
| `S`        | 00        | Matches seconds                                               |
| `s`        | .500      | Matches milliseconds                                          |
| `e`        | Mon, Tues | Matches abbreviated days of the week                          |
| `E`        | Monday    | Matches full name days of the week                            |
| `p`        | AM        | Matches AM/PM (case-insensitive)                              |
| `yyyymmdd` | 19960101  | Matches fixed-width year, month, and day                      |

Characters not listed above are normally treated as delimiters between date and time slots.
For example a `dt` string of "1996-01-15T00:00:00.0" would have a `format` string like
"y-m-dTH:M:S.s". If you need to use a code character as a delimiter you can escape it using
backslash. The date "1995y01m" would have the format "y\\ym\\m".

Note that 12:00AM corresponds 00:00 (midnight), and 12:00PM corresponds to 12:00 (noon).
When parsing a time with a `p` specifier, any hour (either `H` or `I`) is interpreted as
as a 12-hour clock, so the `I` code is mainly useful for output.

Creating a DateFormat object is expensive. Whenever possible, create it once and use it many times
or try the [`dateformat""`](@ref @dateformat_str) string macro. Using this macro creates the DateFormat
object once at macro expansion time and reuses it later. There are also several [pre-defined formatters](@ref
Common-Date-Formatters), listed later.

See [`DateTime`](@ref) and [`format`](@ref) for how to use a DateFormat object to parse and write Date strings
respectively.
"""
function DateFormat(f::AbstractString, locale::DateLocale=ENGLISH)
    tokens = AbstractDateToken[]
    prev = ()
    prev_offset = 1

    # To understand this block, please see the comments attached to the definitions of
    # DATEFORMAT_REGEX_LOCK, DATEFORMAT_REGEX_HASH, and DATEFORMAT_REGEX_CACHE.
    lock(DATEFORMAT_REGEX_LOCK)
    try
        dateformat_regex_hash = hash(keys(CONVERSION_SPECIFIERS))
        if dateformat_regex_hash != DATEFORMAT_REGEX_HASH[]
            DATEFORMAT_REGEX_HASH[] = dateformat_regex_hash
            DATEFORMAT_REGEX_CACHE[] = compute_dateformat_regex(CONVERSION_SPECIFIERS)
        end
    finally
        unlock(DATEFORMAT_REGEX_LOCK)
    end

    for m in eachmatch(DATEFORMAT_REGEX_CACHE[], f)
        tran = replace(f[prev_offset:prevind(f, m.offset)], r"\\(.)" => s"\1")

        if !isempty(prev)
            letter, width = prev
            push!(tokens, DatePart{letter}(width, isempty(tran)))
        end

        if !isempty(tran)
            push!(tokens, Delim(length(tran) == 1 ? first(tran) : tran))
        end

        letter = f[m.offset]
        width = length(m.match)

        prev = (letter, width)
        prev_offset = m.offset + width
    end

    tran = replace(f[prev_offset:lastindex(f)], r"\\(.)" => s"\1")

    if !isempty(prev)
        letter, width = prev
        push!(tokens, DatePart{letter}(width, false))
    end

    if !isempty(tran)
        push!(tokens, Delim(length(tran) == 1 ? first(tran) : tran))
    end

    tokens_tuple = (tokens...,)
    return DateFormat{Symbol(f),typeof(tokens_tuple)}(tokens_tuple, locale)
end

function DateFormat(f::AbstractString, locale::AbstractString)
    DateFormat(f, LOCALES[locale])
end

function Base.show(io::IO, df::DateFormat{S,T}) where {S,T}
    print(io, "dateformat\"", S, '"')
end
Base.Broadcast.broadcastable(x::DateFormat) = Ref(x)

"""
    dateformat"Y-m-d H:M:S"

Create a [`DateFormat`](@ref) object. Similar to `DateFormat("Y-m-d H:M:S")`
but creates the DateFormat object once during macro expansion.

See [`DateFormat`](@ref) for details about format specifiers.
"""
macro dateformat_str(str::String)
    DateFormat(str)
end

# Standard formats

"""
    Dates.ISODateTimeFormat

Describes the ISO8601 formatting for a date and time. This is the default value for `Dates.format`
of a `DateTime`.

# Examples
```jldoctest
julia> Dates.format(DateTime(2018, 8, 8, 12, 0, 43, 1), ISODateTimeFormat)
"2018-08-08T12:00:43.001"
```
"""
const ISODateTimeFormat = DateFormat("yyyy-mm-dd\\THH:MM:SS.s")
default_format(::Type{DateTime}) = ISODateTimeFormat

"""
    Dates.ISODateFormat

Describes the ISO8601 formatting for a date. This is the default value for `Dates.format` of a `Date`.

# Examples
```jldoctest
julia> Dates.format(Date(2018, 8, 8), ISODateFormat)
"2018-08-08"
```
"""
const ISODateFormat = DateFormat("yyyy-mm-dd")
default_format(::Type{Date}) = ISODateFormat

"""
    Dates.ISOTimeFormat

Describes the ISO8601 formatting for a time. This is the default value for `Dates.format` of a `Time`.

# Examples
```jldoctest
julia> Dates.format(Time(12, 0, 43, 1), ISOTimeFormat)
"12:00:43.001"
```
"""
const ISOTimeFormat = DateFormat("HH:MM:SS.s")
default_format(::Type{Time}) = ISOTimeFormat

"""
    Dates.RFC1123Format

Describes the RFC1123 formatting for a date and time.

# Examples
```jldoctest
julia> Dates.format(DateTime(2018, 8, 8, 12, 0, 43, 1), RFC1123Format)
"Wed, 08 Aug 2018 12:00:43"
```
"""
const RFC1123Format = DateFormat("e, dd u yyyy HH:MM:SS")


### API

const Locale = Union{DateLocale, String}

"""
    DateTime(dt::AbstractString, format::AbstractString; locale="english")

Construct a `DateTime` by parsing the `dt` date time string following the
pattern given in the `format` string (see [`DateFormat`](@ref)  for syntax).

!!! note
    This method creates a `DateFormat` object each time it is called. It is recommended
    that you create a [`DateFormat`](@ref) object instead and use that as the second
    argument to avoid performance loss when using the same format repeatedly.

# Examples
```jldoctest
julia> DateTime("2020-01-01", "yyyy-mm-dd")
2020-01-01T00:00:00

julia> a = ("2020-01-01", "2020-01-02");

julia> [DateTime(d, dateformat"yyyy-mm-dd") for d ∈ a] # preferred
2-element Vector{DateTime}:
 2020-01-01T00:00:00
 2020-01-02T00:00:00
```
"""
function DateTime(dt::AbstractString, format::AbstractString; locale::Locale=ENGLISH)
    return parse(DateTime, dt, DateFormat(format, locale))
end

"""
    DateTime(dt::AbstractString, df::DateFormat=ISODateTimeFormat)

Construct a `DateTime` by parsing the `dt` date time string following the
pattern given in the [`DateFormat`](@ref) object, or $ISODateTimeFormat if omitted.

Similar to `DateTime(::AbstractString, ::AbstractString)` but more efficient when
repeatedly parsing similarly formatted date time strings with a pre-created
`DateFormat` object.
"""
DateTime(dt::AbstractString, df::DateFormat=ISODateTimeFormat) = parse(DateTime, dt, df)

"""
    Date(d::AbstractString, format::AbstractString; locale="english")

Construct a `Date` by parsing the `d` date string following the pattern given
in the `format` string (see [`DateFormat`](@ref) for syntax).

!!! note
    This method creates a `DateFormat` object each time it is called. It is recommended
    that you create a [`DateFormat`](@ref) object instead and use that as the second
    argument to avoid performance loss when using the same format repeatedly.

# Examples
```jldoctest
julia> Date("2020-01-01", "yyyy-mm-dd")
2020-01-01

julia> a = ("2020-01-01", "2020-01-02");

julia> [Date(d, dateformat"yyyy-mm-dd") for d ∈ a] # preferred
2-element Vector{Date}:
 2020-01-01
 2020-01-02
```
"""
function Date(d::AbstractString, format::AbstractString; locale::Locale=ENGLISH)
    parse(Date, d, DateFormat(format, locale))
end

"""
    Date(d::AbstractString, df::DateFormat=ISODateFormat)

Construct a `Date` by parsing the `d` date string following the
pattern given in the [`DateFormat`](@ref) object, or $ISODateFormat if omitted.

Similar to `Date(::AbstractString, ::AbstractString)` but more efficient when
repeatedly parsing similarly formatted date strings with a pre-created
`DateFormat` object.
"""
Date(d::AbstractString, df::DateFormat=ISODateFormat) = parse(Date, d, df)

"""
    Time(t::AbstractString, format::AbstractString; locale="english")

Construct a `Time` by parsing the `t` time string following the pattern given
in the `format` string (see [`DateFormat`](@ref) for syntax).

!!! note
    This method creates a `DateFormat` object each time it is called. It is recommended
    that you create a [`DateFormat`](@ref) object instead and use that as the second
    argument to avoid performance loss when using the same format repeatedly.

# Examples
```jldoctest
julia> Time("12:34pm", "HH:MMp")
12:34:00

julia> a = ("12:34pm", "2:34am");

julia> [Time(d, dateformat"HH:MMp") for d ∈ a] # preferred
2-element Vector{Time}:
 12:34:00
 02:34:00
```
"""
function Time(t::AbstractString, format::AbstractString; locale::Locale=ENGLISH)
    parse(Time, t, DateFormat(format, locale))
end

"""
    Time(t::AbstractString, df::DateFormat=ISOTimeFormat)

Construct a `Time` by parsing the `t` date time string following the
pattern given in the [`DateFormat`](@ref) object, or $ISOTimeFormat if omitted.

Similar to `Time(::AbstractString, ::AbstractString)` but more efficient when
repeatedly parsing similarly formatted time strings with a pre-created
`DateFormat` object.
"""
Time(t::AbstractString, df::DateFormat=ISOTimeFormat) = parse(Time, t, df)

@generated function format(io::IO, dt::TimeType, fmt::DateFormat{<:Any,T}) where T
    N = fieldcount(T)
    quote
        ts = fmt.tokens
        loc = fmt.locale
        Base.@nexprs $N i -> format(io, ts[i], dt, loc)
    end
end

function format(dt::TimeType, fmt::DateFormat, bufsize=12)
    # preallocate to reduce resizing
    io = IOBuffer(Vector{UInt8}(undef, bufsize), read=true, write=true)
    format(io, dt, fmt)
    String(io.data[1:io.ptr - 1])
end


"""
    format(dt::TimeType, format::AbstractString; locale="english")::AbstractString

Construct a string by using a `TimeType` object and applying the provided `format`. The
following character codes can be used to construct the `format` string:

| Code       | Examples  | Comment                                                      |
|:-----------|:----------|:-------------------------------------------------------------|
| `y`        | 6         | Numeric year with a fixed width                              |
| `Y`        | 1996      | Numeric year with a minimum width                            |
| `m`        | 1, 12     | Numeric month with a minimum width                           |
| `u`        | Jan       | Month name shortened to 3-chars according to the `locale`    |
| `U`        | January   | Full month name according to the `locale` keyword            |
| `d`        | 1, 31     | Day of the month with a minimum width                        |
| `H`        | 0, 23     | Hour (24-hour clock) with a minimum width                    |
| `M`        | 0, 59     | Minute with a minimum width                                  |
| `S`        | 0, 59     | Second with a minimum width                                  |
| `s`        | 000, 500  | Millisecond with a minimum width of 3                        |
| `e`        | Mon, Tue  | Abbreviated days of the week                                 |
| `E`        | Monday    | Full day of week name                                        |

The number of sequential code characters indicate the width of the code. A format of
`yyyy-mm` specifies that the code `y` should have a width of four while `m` a width of two.
Codes that yield numeric digits have an associated mode: fixed-width or minimum-width.
The fixed-width mode left-pads the value with zeros when it is shorter than the specified
width and truncates the value when longer. Minimum-width mode works the same as fixed-width
except that it does not truncate values longer than the width.

When creating a `format` you can use any non-code characters as a separator. For example to
generate the string "1996-01-15T00:00:00" you could use `format`: "yyyy-mm-ddTHH:MM:SS".
Note that if you need to use a code character as a literal you can use the escape character
backslash. The string "1996y01m" can be produced with the format raw"yyyy\\ymm\\m".
"""
function format(dt::TimeType, f::AbstractString; locale::Locale=ENGLISH)
    format(dt, DateFormat(f, locale))
end

# show
function Base.print(io::IO, dt::DateTime)
    str = if millisecond(dt) == 0
        format(dt, dateformat"YYYY-mm-dd\THH:MM:SS", 19)
    else
        format(dt, dateformat"YYYY-mm-dd\THH:MM:SS.sss", 23)
    end
    print(io, str)
end

function Base.print(io::IO, dt::Date)
    # don't use format - bypassing IOBuffer creation
    # saves a bit of time here.
    y,m,d = yearmonthday(value(dt))
    yy = y < 0 ? @sprintf("%05i", y) : lpad(y, 4, "0")
    mm = lpad(m, 2, "0")
    dd = lpad(d, 2, "0")
    print(io, "$yy-$mm-$dd")
end

for date_type in (:Date, :DateTime)
    # Human readable output (i.e. "2012-01-01")
    @eval Base.show(io::IO, ::MIME"text/plain", dt::$date_type) = print(io, dt)
    # Parsable output (i.e. Date("2012-01-01"))
    @eval Base.show(io::IO, dt::$date_type) = print(io, typeof(dt), "(\"", dt, "\")")
    # Parsable output will have type info displayed, thus it is implied
    @eval Base.typeinfo_implicit(::Type{$date_type}) = true
end
