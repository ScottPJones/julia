# This file is a part of Julia. License is MIT: http://julialang.org/license

include("helpdb/Base.jl")
include("helpdb/statistics.jl")
BUILD_LINALG && include("helpdb/LinAlg.jl")
BUILD_DATES  && include("helpdb/Dates.jl")
BUILD_MMAP   && include("helpdb/Mmap.jl")
