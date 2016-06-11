# This file is a part of Julia. License is MIT: http://julialang.org/license

include("helpdb/Base.jl")
include("helpdb/statistics.jl")
@static if Build.LINALG ; include("helpdb/LinAlg.jl") ; end
@static if Build.DATES ;  include("helpdb/Dates.jl") ; end
@static if Build.MMAP ;   include("helpdb/Mmap.jl") ; end
