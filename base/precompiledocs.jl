# This file is a part of Julia. License is MIT: http://julialang.org/license

# Speed up repl help
sprint(Markdown.term, @doc mean)
sprint(Docs.repl_search, "mean")
sprint(Docs.repl_corrections, "meen")
