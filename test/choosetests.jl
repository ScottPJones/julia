# This file is a part of Julia. License is MIT: http://julialang.org/license

@doc """

`tests, net_on = choosetests(choices)` selects a set of tests to be
run. `choices` should be a vector of test names; if empty or set to
`["all"]`, all tests are selected.

This function also supports "test collections": specifically, "linalg"
 refers to collections of tests in the correspondingly-named
directories.

Upon return, `tests` is a vector of fully-expanded test names, and
`net_on` is true if networking is available (required for some tests).
""" ->
function choosetests(choices = [])
    testnames = [
        "subarray", "core", "inference", "keywordargs", "numbers",
        "printf", "char", "string", "triplequote", "unicode",
        "dict", "hashing", "iobuffer", "staged",
        "arrayops", "tuple", "reduce", "reducedim", "random",
        "abstractarray", "intfuncs", "simdloop", "vecelement",
        "bitarray", "copy", "math", "functional",
        "operators", "path", "ccall", "parse", "loading",
        "sorting", "statistics", "spawn", "backtrace",
        "priorityqueue", "file", "read", "version",
        "pollfd", "broadcast", "socket",
        "floatapprox", "datafmt", "reflection", "regex",
        "combinatorics", "sysinfo", "rounding", "ranges", "mod2pi",
        "euler", "show", "replutil", "sets", "test", "goto", "llvmcall", "grisu",
        "nullable", "meta", "stacktraces", "base64", "serialize", "misc",
        "enums", "cmdlineargs", "i18n", "workspace", "libdl", "int",
        "checked", "intset", "floatfuncs", "compile", "inline",
        "boundscheck", "vecelement", "error", "ambiguous"
    ]

    Base.BUILD_DSP    && push!(testnames, "fft", "dsp")
    Base.BUILD_LINALG && push!(testnames, "linalg", "blas")
    Base.BUILD_SPARSE && push!(testnames, "sparse")
    Base.BUILD_BIGINT && push!(testnames, "bigint", "bignumbers")
    Base.BUILD_BIGFLT && push!(testnames, "mpfr")
    Base.BUILD_DATES  && push!(testnames, "dates")
    Base.BUILD_PKG    && push!(testnames, "libgit2", "resolve")
    Base.BUILD_REPL   && push!(testnames, "lineedit", "replcompletions", "repl")
    Base.BUILD_MMAP   && push!(testnames, "mmap")
    Base.BUILD_DOCS   && push!(testnames, "docs", "markdown")
    Base.BUILD_FULL   && push!(testnames, "fastmath")
    Base.BUILD_FLOAT16  && push!(testnames, "float16")
    Base.BUILD_PROFILER && push!(testnames, "profile")
    Base.BUILD_COMPLEX  && push!(testnames, "complex")
    Base.BUILD_THREADS  && push!(testnames, "threads")
    Base.BUILD_PARALLEL && push!(testnames, "parallel")

    if isdir(joinpath(JULIA_HOME, Base.DOCDIR, "examples"))
        push!(testnames, "examples")
    end

    tests = []
    skip_tests = []

    for (i, t) in enumerate(choices)
        if t == "--skip"
            skip_tests = choices[i + 1:end]
            break
        else
            push!(tests, t)
        end
    end

    if tests == ["all"] || isempty(tests)
        tests = testnames
    end

    linalgtests = ["linalg/triangular", "linalg/qr", "linalg/dense",
                   "linalg/matmul", "linalg/schur", "linalg/special",
                   "linalg/eigen", "linalg/bunchkaufman", "linalg/svd",
                   "linalg/lapack", "linalg/tridiag", "linalg/bidiag",
                   "linalg/diagonal", "linalg/pinv", "linalg/givens",
                   "linalg/cholesky", "linalg/lu", "linalg/symmetric",
                   "linalg/generic", "linalg/uniformscaling", "linalg/lq",
                   "linalg/hessenberg"]
    if Base.USE_GPL_LIBS
        push!(linalgtests, "linalg/arnoldi")
    end

    if "linalg" in skip_tests
        filter!(x -> (x != "linalg" && !(x in linalgtests)), tests)
    elseif "linalg" in tests
        # specifically selected case
        filter!(x -> x != "linalg", tests)
        prepend!(tests, linalgtests)
    end

    net_required_for = ["socket", "parallel"]
    net_on = true
    try
        getipaddr()
    catch
        warn("Networking unavailable: Skipping tests [" * join(net_required_for, ", ") * "]")
        net_on = false
    end

    if ccall(:jl_running_on_valgrind,Cint,()) != 0 && "rounding" in tests
        warn("Running under valgrind: Skipping rounding tests")
        filter!(x -> x != "rounding", tests)
    end

    if !net_on
        filter!(x -> !(x in net_required_for), tests)
    end

    filter!(x -> !(x in skip_tests), tests)

    tests, net_on
end
