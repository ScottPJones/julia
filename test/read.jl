# This file is a part of Julia. License is MIT: https://julialang.org/license

using Random, Sockets

mktempdir() do dir

tasks = []

# Create test file
filename = joinpath(dir, "file.txt")
text = "C1,C2\n1,2\na,b\n"

# List of IO producers
l = Vector{Tuple{AbstractString,Function}}()


# File
io = (text) -> begin
    write(filename, text)
    Base.Filesystem.open(filename, Base.Filesystem.JL_O_RDONLY)
end
s = io(text)
@test isa(s, IO)
@test isa(s, Base.Filesystem.File)
close(s)
push!(l, ("File", io))


# IOStream
io = (text) -> begin
    write(filename, text)
    open(filename)
end
s = io(text)
@test isa(s, IO)
@test isa(s, IOStream)
close(s)
push!(l, ("IOStream", io))


# IOBuffer
io = (text)->IOBuffer(text)
s = io(text)
@test isa(s, IO)
@test isa(s, IOBuffer)
close(s)
push!(l, ("IOBuffer", io))


function run_test_server(srv, text)
    push!(tasks, @async begin
        try
            sock = accept(srv)
            try
                write(sock, text)
            catch e
                if !(isa(e, Base.IOError) && e.code == Base.UV_EPIPE)
                    if !(isa(e, Base.IOError) && e.code == Base.UV_ECONNRESET)
                        rethrow()
                    end
                end
            finally
                close(sock)
            end
        finally
            close(srv)
        end
    end)
    yield()
end


# TCPSocket
io = (text) -> begin
    port, srv = listenany(rand(2000:4000))
    run_test_server(srv, text)
    connect(port)
end
s = io(text)
@test isa(s, IO)
@test isa(s, TCPSocket)
close(s)
push!(l, ("TCPSocket", io))


# PipeEndpoint
io = (text) -> begin
    a = "\\\\.\\pipe\\uv-test-$(randstring(6))"
    b = joinpath(dir, "socket-$(randstring(6))")
    socketname = Sys.iswindows() ? a : b
    srv = listen(socketname)
    run_test_server(srv, text)
    connect(socketname)
end
s = io(text)
@test isa(s, IO)
@test isa(s, Base.PipeEndpoint)
close(s)
push!(l, ("PipeEndpoint", io))


# Pipe (#14747)
io = (text) -> begin
    write(filename, text)
    # we can skip using shell_escape_wincmd, since ", ^, and % aren't legal in
    # a filename, so unconditionally wrapping in " is sufficient (okay, that's
    # a lie, since ^ and % actually are legal, but DOS is broken)
    if Sys.iswindows()
        cmd = Cmd(["cmd.exe", "/c type \"$(replace(filename, '/' => '\\'))\""])
        cmd = Cmd(cmd; windows_verbatim=true)
        cmd = pipeline(cmd, stderr=devnull)
    else
        cmd = `cat $filename`
    end
    open(cmd)
end
s = io(text)
@test isa(s, IO)
@test isa(s, Base.Process)
close(s)
push!(l, ("Process", io))


open_streams = []
function cleanup()
    for s_ in open_streams
        close(s_)
    end
    empty!(open_streams)
    for tsk in tasks
        wait(tsk)
    end
    empty!(tasks)
end


verbose = false

for (name, f) in l
    local function io(text=text)
        local s = f(text)
        push!(open_streams, s)
        return s
    end

    verbose && println("$name readuntil...")
    for (t, s, m, kept) in [
            ("a", "", "", ""),
            ("a", "ab", "a", "a"),
            ("b", "ab", "b", "b"),
            ("α", "αγ", "α", "α"),
            ("ab", "abc", "ab", "ab"),
            ("bc", "abc", "bc", "bc"),
            ("αβ", "αβγ", "αβ", "αβ"),
            ("aaabc", "ab", "aa", "aaab"),
            ("aaabc", "b", "aaa", "aaab"),
            ("aaabc", "ac", "aaabc", "aaabc"),
            ("aaabc", "aab", "a", "aaab"),
            ("aaabc", "aac", "aaabc", "aaabc"),
            ("αααβγ", "αβ", "αα", "αααβ"),
            ("αααβγ", "β", "ααα", "αααβ"),
            ("αααβγ", "ααβ", "α", "αααβ"),
            ("αααβγ", "αγ", "αααβγ", "αααβγ"),
            ("barbarbarians", "barbarian", "bar", "barbarbarian"),
            ("abcaabcaabcxl", "abcaabcx", "abca", "abcaabcaabcx"),
            ("abbaabbaabbabbaax", "abbaabbabbaax", "abba", "abbaabbaabbabbaax"),
            ("abbaabbabbaabbaabbabbaax", "abbaabbabbaax", "abbaabbabba", "abbaabbabbaabbaabbabbaax"),
            ('a'^500 * 'x' * "bbbb", "x", 'a'^500, 'a'^500 * 'x'),
           ]
        local t, s, m, kept
        @test readuntil(io(t), s) == m
        @test readuntil(io(t), s, keep=true) == kept
        if isone(length(s))
            @test readuntil(io(t), first(s)) == m
            @test readuntil(io(t), first(s), keep=true) == kept
        end
        @test readuntil(io(t), SubString(s, firstindex(s))) == m
        @test readuntil(io(t), SubString(s, firstindex(s)), keep=true) == kept
        @test readuntil(io(t), GenericString(s)) == m
        @test readuntil(io(t), GenericString(s), keep=true) == kept
        @test readuntil(io(t), unsafe_wrap(Vector{UInt8},s)) == unsafe_wrap(Vector{UInt8},m)
        @test readuntil(io(t), unsafe_wrap(Vector{UInt8},s), keep=true) == unsafe_wrap(Vector{UInt8},kept)
        @test readuntil(io(t), collect(s)::Vector{Char}) == Vector{Char}(m)
        @test readuntil(io(t), collect(s)::Vector{Char}, keep=true) == Vector{Char}(kept)

        buf = IOBuffer()
        @test String(take!(copyuntil(buf, io(t), s))) == m
        @test String(take!(copyuntil(buf, io(t), s, keep=true))) == kept
        file = tempname()
        for (k,m) in ((false, m), (true, kept))
            open(file, "w") do f
                @test f == copyuntil(f, io(t), s, keep=k)
            end
            @test read(file, String) == m
        end
        rm(file)
    end
    cleanup()

    write(filename, text)

    verbose && println("$name read...")
    @test read(io(), UInt8) == read(IOBuffer(text), UInt8)
    @test read(io(), UInt8) == read(filename, UInt8)
    @test read(io(), Int) == read(IOBuffer(text), Int)
    @test read(io(), Int) == read(filename,Int)
    s1 = io()
    s2 = IOBuffer(text)
    @test read!(s1, Vector{UInt32}(undef, 2)) == read!(s2, Vector{UInt32}(undef, 2))
    @test !eof(s1)
    @test read!(s1, Vector{UInt8}(undef, 5)) == read!(s2, Vector{UInt8}(undef, 5))
    @test !eof(s1)
    @test read!(s1, Vector{UInt8}(undef, 1)) == read!(s2, Vector{UInt8}(undef, 1))
    @test eof(s1)
    @test_throws EOFError read(s1, UInt8)
    @test eof(s1)
    close(s1)
    close(s2)

    verbose && println("$name eof...")
    n = length(text) - 1
    @test read!(io(), Vector{UInt8}(undef, n)) ==
          read!(IOBuffer(text), Vector{UInt8}(undef, n))
    @test (s = io(); read!(s, Vector{UInt8}(undef, n)); !eof(s))
    n = length(text)
    @test read!(io(), Vector{UInt8}(undef, n)) ==
          read!(IOBuffer(text), Vector{UInt8}(undef, n))
    @test (s = io(); read!(s, Vector{UInt8}(undef, n)); eof(s))
    n = length(text) + 1
    @test_throws EOFError read!(io(), Vector{UInt8}(undef, n))
    @test_throws EOFError read!(io(), Vector{UInt8}(undef, n))

    old_text = text
    cleanup()

    for text_ in [
        old_text,
        String(Char['A' + i % 52 for i in 1:(div(Base.SZ_UNBUFFERED_IO,2))]),
        String(Char['A' + i % 52 for i in 1:(    Base.SZ_UNBUFFERED_IO -1)]),
        String(Char['A' + i % 52 for i in 1:(    Base.SZ_UNBUFFERED_IO   )]),
        String(Char['A' + i % 52 for i in 1:(    Base.SZ_UNBUFFERED_IO +1)])
    ]
        text = text_
        write(filename, text)

        verbose && println("$name read(io, String)...")
        @test read(io(), String) == text

        @test read(io(), String) == read(filename, String)


        verbose && println("$name read...")
        @test read(io()) == unsafe_wrap(Vector{UInt8},text)

        @test read(io()) == read(filename)

        cleanup()


        verbose && println("$name readbytes!...")
        l = length(text)
        for n = [1, 2, l-2, l-1, l, l+1, l+2]
            a1 = Vector{UInt8}(undef, n)
            a2 = Vector{UInt8}(undef, n)
            s1 = io()
            s2 = IOBuffer(text)
            n1 = readbytes!(s1, a1)
            n2 = readbytes!(s2, a2)
            @test n1 == n2
            @test length(a1) == length(a2)
            let l = min(l, n)
                @test a1[1:l] == a2[1:l]
            end
            @test n <= length(text) || eof(s1)
            @test n <= length(text) || eof(s2)

            cleanup()
        end

        # Test growing output array
        let x = UInt8[],
            io = io()
            n = readbytes!(io, x)
            @test n == 0
            @test isempty(x)
            n = readbytes!(io, x, typemax(Int))
            @test n == length(x)
            @test x == codeunits(text)
            cleanup()
        end

        verbose && println("$name read!...")
        l = length(text)
        for n = [1, 2, l-2, l-1, l]
            @test read!(io(), Vector{UInt8}(undef, n)) ==
                  read!(IOBuffer(text), Vector{UInt8}(undef, n))
            @test read!(io(), Vector{UInt8}(undef, n)) ==
                  read!(filename, Vector{UInt8}(undef, n))

            cleanup()
        end
        @test_throws EOFError read!(io(), Vector{UInt8}(undef, length(text)+1))


        verbose && println("$name readuntil...")
        for keep in [false, true]
            @test readuntil(io(), '\n', keep=keep) == readuntil(IOBuffer(text),'\n', keep=keep)
            @test readuntil(io(), '\n', keep=keep) == readuntil(filename,'\n', keep=keep)
            @test readuntil(io(), "\n", keep=keep) == readuntil(IOBuffer(text),"\n", keep=keep)
            @test readuntil(io(), "\n", keep=keep) == readuntil(filename,"\n", keep=keep)
            @test readuntil(io(), ',', keep=keep)  == readuntil(IOBuffer(text),',', keep=keep)
            @test readuntil(io(), ',', keep=keep)  == readuntil(filename,',', keep=keep)
        end

        cleanup()

        verbose && println("$name readline...")
        file = tempname()
        for lineending in ("\n", "\r\n", "")
            kept = "foo bar" * lineending
            t = isempty(lineending) ? "foo bar" : kept * "baz\n"
            write(file, t)
            @test readline(io(t)) == readline(file) == "foo bar"
            @test readline(io(t), keep=true) == readline(file, keep=true) == kept

            @test String(take!(copyline(IOBuffer(), file))) == "foo bar"
            @test String(take!(copyline(IOBuffer(), file, keep=true))) == kept

            cleanup()

            buf = IOBuffer()
            @test buf === copyline(buf, io(t))
            @test String(take!(buf)) == "foo bar"
            @test String(take!(copyline(buf, file, keep=true))) == kept
            for keep in (true, false)
                open(file, "w") do f
                    @test f === copyline(f, io(t), keep=keep)
                end
                @test read(file, String) == (keep ? kept : "foo bar")
            end

            cleanup()

            write(file, lineending)
            @test readline(IOBuffer(lineending)) == ""
            @test readline(IOBuffer(lineending), keep=true) == lineending
            @test String(take!(copyline(IOBuffer(), IOBuffer(lineending)))) == ""
            @test String(take!(copyline(IOBuffer(), IOBuffer(lineending), keep=true))) == lineending
            @test readline(file) == ""
            @test readline(file, keep=true) == lineending
            @test String(take!(copyline(IOBuffer(), file))) == ""
            @test String(take!(copyline(IOBuffer(), file, keep=true))) == lineending

            cleanup()
        end
        rm(file)

        verbose && println("$name readlines...")
        @test readlines(io(), keep=true) == readlines(IOBuffer(text), keep=true)
        @test readlines(io(), keep=true) == readlines(filename, keep=true)
        @test readlines(io()) == readlines(IOBuffer(text))
        @test readlines(io()) == readlines(filename)
        @test collect(eachline(io(), keep=true)) == collect(eachline(IOBuffer(text), keep=true))
        @test collect(eachline(io(), keep=true)) == collect(eachline(filename, keep=true))
        @test collect(eachline(io())) == collect(eachline(IOBuffer(text)))
        @test collect(@inferred(eachline(io()))) == collect(@inferred(eachline(filename))) #20351
        if try; seekend(io()); true; catch; false; end # reverse iteration only supports seekable streams
            for keep in (true, false)
                lines = readlines(io(); keep)
                @test last(lines) == last(eachline(io(); keep))
                @test last(lines,2) == last(eachline(io(); keep),2)
                @test reverse!(lines) == collect(Iterators.reverse(eachline(io(); keep))) == collect(Iterators.reverse(eachline(IOBuffer(text); keep)))
            end
        end

        cleanup()

        verbose && println("$name readeach...")
        @test collect(readeach(io(), Char)) == Vector{Char}(text)
        @test collect(readeach(io(), UInt8)) == Vector{UInt8}(text)

        cleanup()

        verbose && println("$name countlines...")
        @test countlines(io()) == countlines(IOBuffer(text))
    end

    text = old_text
    write(filename, text)

    if !isa(io(), Union{Base.PipeEndpoint, Base.AbstractPipe, TCPSocket})
        verbose && println("$name position...")
        @test (s = io(); read!(s, Vector{UInt8}(undef, 4)); position(s)) == 4

        verbose && println("$name seek...")
        for n = 0:length(text)-1
            @test readlines(seek(io(), n)) == readlines(seek(IOBuffer(text), n))
            cleanup()
        end
        verbose && println("$name skip...")
        for n = 0:length(text)-1
            @test readlines(seek(io(), n)) == readlines(seek(IOBuffer(text), n))
            @test readlines(skip(io(), n)) == readlines(skip(IOBuffer(text), n))
            cleanup()
        end
        verbose && println("$name seekend...")
        @test read(seekend(io()), String) == ""
    end


    verbose && println("$name write(::IOStream, ...)")
    to = open("$filename.to", "w")
    write(to, io())
    close(to)
    @test read("$filename.to", String) == text

    verbose && println("$name write(filename, ...)")
    write("$filename.to", io())
    @test read("$filename.to", String) == text

    verbose && println("$name write(::IOBuffer, ...)")
    to = IOBuffer(Vector{UInt8}(codeunits(text)), read=false, write=true)
    write(to, io())
    @test String(take!(to)) == text

    cleanup()
end

function test_read_nbyte()
    fn = tempname()
    # Write one byte. One byte read should work once
    # but 2-byte read should throw EOFError.
    open(fn, "w+") do f
        write(f, 0x55)
        flush(f)
        seek(f, 0)
        @test read(f, UInt8) == 0x55
        @test_throws EOFError read(f, UInt8)
        seek(f, 0)
        @test_throws EOFError read(f, UInt16)
    end
    # Write 2 more bytes. Now 2-byte read should work once
    # but 4-byte read should fail with EOFError.
    open(fn, "a+") do f
        write(f, 0x4444)
        flush(f)
        seek(f, 0)
        @test read(f, UInt16) == 0x4455
        @test_throws EOFError read(f, UInt16)
        seek(f, 0)
        @test_throws EOFError read(f, UInt32)
    end
    # Write 4 more bytes. Now 4-byte read should work once
    # but 8-byte read should fail with EOFError.
    open(fn, "a+") do f
        write(f, 0x33333333)
        flush(f)
        seek(f, 0)
        @test read(f, UInt32) == 0x33444455
        @test_throws EOFError read(f, UInt32)
        seek(f, 0)
        @test_throws EOFError read(f, UInt64)
    end
    # Writing one more byte should allow an 8-byte
    # read to proceed.
    open(fn, "a+") do f
        write(f, 0x22)
        flush(f)
        seek(f, 0)
        @test read(f, UInt64) == 0x2233333333444455
    end
    rm(fn)
end
test_read_nbyte()


# devnull
@test !isreadable(devnull)
@test iswritable(devnull)
@test isopen(devnull)
@test write(devnull, 0xff) === 1
@test write(devnull, Int32(1234)) === 4
@test_throws EOFError read(devnull, UInt8)
@test close(devnull) === nothing
@test flush(devnull) === nothing
@test eof(devnull)
@test print(devnull, "go to /dev/null") === nothing


let s = "qwerty"
    @test read(IOBuffer(s)) == codeunits(s)
    @test read(IOBuffer(s), 10) == codeunits(s)
    @test read(IOBuffer(s), 1) == codeunits(s)[1:1]
end


# Filesystem.File
f = joinpath(dir, "test.txt")
open(io->write(io, "123"), f, "w")
f1 = open(f)
f2 = Base.Filesystem.open(f, Base.Filesystem.JL_O_RDONLY)
@test read(f1, UInt8) == read(f2, UInt8) == UInt8('1')
@test read(f1, UInt8) == read(f2, UInt8) == UInt8('2')
@test read(f1, UInt8) == read(f2, UInt8) == UInt8('3')
@test_throws EOFError read(f1, UInt8)
@test_throws EOFError read(f2, UInt8)
close(f1)
close(f2)

a = UInt8[0,0,0]
f1 = open(f)
f2 = Base.Filesystem.open(f, Base.Filesystem.JL_O_RDONLY)
@test read!(f1, a) == read!(f2, a) == UInt8['1','2','3']
@test_throws EOFError read!(f1, a)
@test_throws EOFError read!(f2, a)
close(f1)
close(f2)

a = UInt8[0,0,0,0]
f1 = open(f)
f2 = Base.Filesystem.open(f, Base.Filesystem.JL_O_RDONLY)
@test_throws EOFError read!(f1, a)
@test_throws EOFError read!(f2, a)
close(f1)
close(f2)
rm(f)

io = Base.Filesystem.open(f, Base.Filesystem.JL_O_WRONLY | Base.Filesystem.JL_O_CREAT | Base.Filesystem.JL_O_EXCL, 0o000)
@test write(io, "abc") == 3
close(io)
if !Sys.iswindows() && Libc.geteuid() != 0 # root user
    # msvcrt _wchmod documentation states that all files are readable,
    # so we don't test that it correctly set the umask on windows
    @test_throws SystemError open(f)
    @test_throws Base.IOError Base.Filesystem.open(f, Base.Filesystem.JL_O_RDONLY)
else
    Sys.iswindows() || @warn "File permissions tests skipped due to running tests as root (not recommended)"
    close(open(f))
end
chmod(f, 0o400)
f1 = open(f)
f2 = Base.Filesystem.open(f, Base.Filesystem.JL_O_RDONLY)
for i = 1:2
    @test !eof(f1)
    @test !eof(f2)
    @test position(f1) == 0
    @test position(f2) == 0
    @test read(f1, String) == read(f2, String) == "abc"
    @test read(f1, String) == read(f2, String) == ""
    @test position(f1) == 3
    @test position(f2) == 3
    @test eof(f1)
    @test eof(f2)
    @test seekstart(f1) == f1
    @test seekstart(f2) == f2
end
@test seekend(f1) == f1
@test seekend(f2) == f2
@test eof(f1)
@test eof(f2)
@test skip(f1, -2) == f1
@test skip(f2, -2) == f2
@test position(f1) == 1
@test position(f2) == 1
@test_throws SystemError skip(f1, -2)
@test_throws SystemError skip(f2, -2)
@test position(f1) == 1
@test position(f2) == 1
@test skip(f1, 300) == f1
@test skip(f2, 300) == f2
@test position(f1) == 301
@test position(f2) == 301
@test eof(f1)
@test eof(f2)
@test_throws ArgumentError write(f1, '*')
@test_throws Base.IOError write(f2, '*')
close(f1)
close(f2)
@test eof(f1)
@test_throws Base.IOError eof(f2)
if Libc.geteuid() != 0 # root user
    @test_throws SystemError open(f, "r+")
    @test_throws Base.IOError Base.Filesystem.open(f, Base.Filesystem.JL_O_RDWR)
else
    @warn "File permissions tests skipped due to running tests as root (not recommended)"
end
chmod(f, 0o600)
f1 = open(f, "r+")
f2 = Base.Filesystem.open(f, Base.Filesystem.JL_O_RDWR)
@test skip(f1, 10) == f1
@test skip(f2, 10) == f2
@test eof(f1)
@test eof(f2)
@test write(f1, '*') == 1
@test flush(f1) === nothing
@test !eof(f2)
@test skip(f2, 1) == f2
@test write(f2, '*') == 1
@test !eof(f1)
@test seekstart(f1) == f1
@test seekstart(f2) == f2
@test read(f1, String) == read(f2, String) == "abc\0\0\0\0\0\0\0**"
close(f1)
close(f2)
rm(f)

end # mktempdir() do dir

@testset "countlines" begin
    @test countlines(IOBuffer("")) == 0
    @test countlines(IOBuffer("\n")) == 1
    @test countlines(IOBuffer("\n"), eol = '\r') == 1
    @test countlines(IOBuffer("\r\r\n\r"), eol = '\r') == 3
    @test countlines(IOBuffer("\n\n\n\n\n\n\n\n\n\n")) == 10
    @test countlines(IOBuffer("\n \n \n \n \n \n \n \n \n \n")) == 10
    @test countlines(IOBuffer("\r\n \r\n \r\n \r\n \r\n")) == 5
    @test countlines(IOBuffer("foo\nbar")) == length(readlines(IOBuffer("foo\nbar"))) == 2
    file = tempname()
    write(file,"Spiffy header\nspectacular first row\neven better 2nd row\nalmost done\n")
    @test countlines(file) == 4
    @test countlines(file, eol = '\r') == 1
    @test countlines(file, eol = '\n') == 4
    rm(file)
end

let p = Pipe()
    Base.link_pipe!(p, reader_supports_async=true, writer_supports_async=true)
    t = @async read(p)
    @sync begin
        @async write(p, zeros(UInt16, 660_000))
        yield() # TODO: need to add an Event to the previous line
        order::UInt16 = 0
        for i = 1:typemax(UInt16)
            @async (order += 1; write(p, order); nothing)
        end
        yield() # TODO: need to add an Event to the previous line
        @async close(p.in)
    end
    s = reinterpret(UInt16, fetch(t))
    @test length(s) == 660_000 + typemax(UInt16)
    @test s[(end - typemax(UInt16)):end] == UInt16.(0:typemax(UInt16))
end

# issue #26419
@test Base.return_types(read, (String, Type{String})) == Any[String]

@testset "read! to view" begin
    x = rand(4, 4)
    y = rand(10)
    z = 1:10
    v = [1.0, 2.0, 3.0, 4.0]
    io = IOBuffer()
    write(io, v)
    flush(io)
    seekstart(io)
    read!(io, @view x[:, 3])
    @test x[:, 3] == v
    x = rand(3, 3)
    seekstart(io)
    read!(io, @view x[1:2, 2:3])
    @test x[1:2, 2:3][:] == v[:]
    seekstart(io)
    read!(io, @view y[4:7])
    @test y[4:7] == v
    seekstart(io)
    @test_throws Base.CanonicalIndexError read!(io, @view z[4:6])
end

# Bulk read from pipe
let p = Pipe()
    data = rand(UInt8, Base.SZ_UNBUFFERED_IO + 100)
    Base.link_pipe!(p, reader_supports_async=true, writer_supports_async=true)
    t = @async write(p.in, data)
    @test read(p.out, UInt8) == data[1]
    data_read = Vector{UInt8}(undef, 10*Base.SZ_UNBUFFERED_IO)
    nread = readbytes!(p.out, data_read, Base.SZ_UNBUFFERED_IO + 50)
    @test nread == Base.SZ_UNBUFFERED_IO + 50
    @test data_read[1:nread] == data[2:nread+1]
    @test read(p.out, 49) == data[end-48:end]
    wait(t)

    closewrite(p)
    @test !isopen(p.in)
    @test isopen(p.out)
    close(p)
    @test !isopen(p.out)
end

@testset "issue #27412" for itr in [eachline(IOBuffer("a")), readeach(IOBuffer("a"), Char)]
    @test !isempty(itr)
    # check that the earlier isempty did not consume the iterator
    @test !isempty(itr)
    first(itr) # consume the iterator
    @test  isempty(itr) # now it is empty
end

@testset "readuntil/copyuntil fallbacks" begin
    # test fallback for generic delim::T
    buf = IOBuffer()
    fib = [1,1,2,3,5,8,13,21]
    write(buf, fib)
    @test readuntil(seekstart(buf), 21) == fib[1:end-1]
    @test readuntil(buf, 21) == Int[]
    @test readuntil(seekstart(buf), 21; keep=true) == fib
    out = IOBuffer()
    @test copyuntil(out, seekstart(buf), 21) === out
    @test reinterpret(Int, take!(out)) == fib[1:end-1]
    @test copyuntil(out, seekstart(buf), 21; keep=true) === out
    @test reinterpret(Int, take!(out)) == fib
end

# more tests for reverse(eachline)
@testset "reverse(eachline)" begin
    lines = vcat(repr.(1:4), ' '^50000 .* repr.(5:10), repr.(11:10^5))
    for lines in (lines, reverse(lines)), finalnewline in (true, false), eol in ("\n", "\r\n")
        buf = IOBuffer(join(lines, eol) * (finalnewline ? eol : ""))
        @test reverse!(collect(Iterators.reverse(eachline(seekstart(buf))))) == lines
        @test last(eachline(seekstart(buf))) == last(lines)
        @test last(eachline(seekstart(buf)),10^4) == last(lines,10^4)
        @test last(eachline(seekstart(buf)),length(lines)*2) == lines
        @test reverse!(collect(Iterators.reverse(eachline(seek(buf, sum(sizeof, lines[1:100]) + 100*sizeof(eol)))))) == lines[101:end]
        @test isempty(Iterators.reverse(eachline(buf)))
    end

    let rempty = Iterators.reverse(eachline(IOBuffer()))
        @test isempty(rempty)
        @test isempty(collect(rempty))
    end

    let buf = IOBuffer("foo\nbar")
        @test readline(buf) == "foo"
        r = Iterators.reverse(eachline(buf))
        line, state = iterate(r)
        @test line == "bar"
        @test Base.isdone(r, state)
        @test Base.isdone(r)
        @test isempty(r) && isempty(collect(r))
    end
end

@testset "Ref API" begin
    io = PipeBuffer()
    @test write(io, Ref{Any}(0xabcd_1234)) === 4
    @test read(io, UInt32) === 0xabcd_1234
    @test_throws ErrorException("write cannot copy from a Ptr") invoke(write, Tuple{typeof(io), Ref{Cvoid}}, io, C_NULL)
    @test_throws ErrorException("write cannot copy from a Ptr") invoke(write, Tuple{typeof(io), Ref{Int}}, io, Ptr{Int}(0))
    @test_throws ErrorException("write cannot copy from a Ptr") invoke(write, Tuple{typeof(io), Ref{Any}}, io, Ptr{Any}(0))
    @test_throws ErrorException("read! cannot copy into a Ptr") read!(io, C_NULL)
    @test_throws ErrorException("read! cannot copy into a Ptr") read!(io, Ptr{Int}(0))
    @test_throws ErrorException("read! cannot copy into a Ptr") read!(io, Ptr{Any}(0))
    @test eof(io)
    @test write(io, C_NULL) === sizeof(Int)
    @test write(io, Ptr{Int}(4)) === sizeof(Int)
    @test write(io, Ptr{Any}(5)) === sizeof(Int)
    @test read!(io, Int[1, 2, 3]) == [0, 4, 5]
    @test eof(io)
end
