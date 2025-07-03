"""
    $(@__MODULE__)

A combinatorics package for Julia based on SmallCollections.jl.

In the examples we assume that the package is loaded together with SmallCollections:
```julia
julia> using SmallCollections, $(@__MODULE__)
```
"""
module SmallCombinatorics

using Base: Fix1, Fix2, Generator

using SmallCollections
using SmallCollections: bitsize, padtail, unsafe_shl, unsafe_lshr,
    blsi, blsr, blsmsk, pdep,
    AbstractFixedOrSmallVector, AbstractFixedOrSmallOrPackedVector

_SmallBitSet(mask::U) where U <: Unsigned = convert(SmallBitSet{U}, mask)

include("partitions.jl")
include("compositions.jl")
include("subsets.jl")
include("permutations.jl")

end # module
