using Documenter
using SmallCollections, SmallCombinatorics

DocMeta.setdocmeta!(SmallCombinatorics, :DocTestSetup, quote
        using SmallCollections, SmallCombinatorics
        # for jldoctest in docstrings
    end; recursive = true)

makedocs(sitename = "SmallCombinatorics.jl",
    modules = [
        SmallCombinatorics,
    ],
    pages = [
        "index.md",
    ],
    format = Documenter.HTML(),
    warnonly = true)
