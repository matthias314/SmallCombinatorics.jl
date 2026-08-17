
export ldm, NoSubsets

using DataStructures: BinaryMaxHeap, AbstractHeap
using StaticArrays: SVector, MVector

struct NoSubsets <: AbstractSet{Any} end
NoSubsets(x...; kw...) = NoSubsets()
Base.push!(ns::NoSubsets, args...) = ns
Base.length(ns::NoSubsets) = 0
Base.union(ns::NoSubsets, args...) = ns
Base.union!(ns::NoSubsets, args...) = ns

function settype(::Val{n}) where n
    n <= 512 && return smallbitsettype(Val(n))
    return BitSet
end

struct LDMNode{T,K,SD}
    value::T
    subsets::NTuple{K,SD}
    subsetsums::NTuple{K,T}
end

value(ld) = ld.value
subs(ld) = ld.subsets
subs(ld, i) = ld.subsets[i]
sums(ld) = ld.subsetsums
sums(ld, i) = ld.subsetsums[i]
Base.isless(a::LDMNode, b::LDMNode) = value(a) < value(b)

function LDMNode{K, SD}(n, idx=nothing) where {K,SD}
    singleton = SD<:Union{SmallBitSet,SmallSet} ?
        idx === nothing ? SD(n) : SD(idx) :
        idx === nothing ? push!(SD(), n) : push!(SD(), idx)

    subsets = ntuple(Val(K)) do i
        i==1 ? singleton : SD()
    end

    sums = ntuple(Val(K)) do i
        i==1 ? n : zero(n)
    end

    LDMNode(n, subsets, sums)
end



"""ldm(numbers, ::Val{K}=Val(2), ::Type=<Bitset type>)

Compute partitions of an indexable set of non-negative numbers, i.e. an `AbstractVector`.
That is, the set of numbers is partitioned into subsets such that
their sums are close.

The Largest Differencing Method (Karmarkar-Karp) is used.

The second argument is the number of partitions expressed as a `Val`, e.g. `Val(2)`
for splitting a set in two disjoint subsets with close sums.

The last argument must be a set type suitable for holding subsets of indices into
`numbers`.  It defaults to
`SmallBitSet{UInt8/16/32/64/128/256/512}` or `BitSet`. If an
unsigned type `T` is supplied, `SmallBitSet{T}` will be used. `NoSubsets` can also be
specified, then the subsets will not be computed, only their sums.
This can save considerable time.

Note that the default choice is not type stable unless `numbers` is a fixed size
vector type like `FixedVector` or `StaticVector`. That is, a suitable set type is chosen
depending on the cardinality of `numbers`. Crossing from data into types
is the definition of type instability. In critical code sections you
should therefore supply a large enough type, or `BitSet`.
That is, use `UInt64` if you know there are no more than 64 numbers in `numbers`.

A struct is returned with 3 fields:
* value: the difference between the largest and smallest subset sum
* subsets: a tuple with sets of indices for each of the subsets
* subsetsums: a tuple with the sums of the subsets


# Examples
```julia
julia> s = [0.3, 0.2, 0.6, 0.23];
julia> ldm(s, Val(2)).subsetsums
(0.73, 0.6)

julia> ldm(1:11).subsets
(SmallBitSet{UInt16}([2, 4, 7, 9, 11]), SmallBitSet{UInt16}([1, 3, 5, 6, 8, 10]))

julia> s = rand(102:3:137, 14)
julia> subsets = ldm(s, Val(3), UInt16).subsets
julia> [[s[i] for i in ind] for ind in subsets]
3-element Vector{Vector{Int64}}:
 [102, 111, 132, 129, 117]
 [135, 105, 102, 120, 129]
 [129, 135, 117, 123]
```
"""
function ldm(numbers::AbstractVector, ::Val{K}=Val(2),
             ::Type{ST}=settype(Val(length(numbers)))) where {K,ST<:AbstractSet}
    Base.require_one_based_indexing(numbers)
    heap = BinaryMaxHeap(LDMNode{K,ST}.(numbers, eachindex(numbers)))
    ldm(heap, Val(K))
end

function ldm(numbers::AbstractSet, ::Val{K}=Val(2),
             ::Type{ST} = typeof(numbers)) where {K,ST<:AbstractSet}
    heap = BinaryMaxHeap(LDMNode{K,ST}.(numbers))
    ldm(heap, Val(K))
end


#=
"""ldm(heap::AbstractHeap{VT}, ::Val{K}) where {K,VT}

Variant of `ldm` which does not take a `numbers` vector, but a (max) heap
with elements of type `VT`.  There must be extensions of the functions
`value(::VT)`, `subs(::VT)`, `subs(::VT, i)` and `sums(::VT)`,
`sums(::VT, i)`, and `Base.isless(::VT, ::VT)`.  And a constructor
`VT(val, subsets, subsetsums)`.

On input the `heap` must contain one entry for each number, with a
singleton subset.

The heap will be modified and emptied (`pop!` and `push!`),
and an element of it will be returned.
This somewhat awkward interface can be of use if one has a heap which can be
reused without allocations.
"""
=#
function ldm(heap::AbstractHeap{VT}, ::Val{K}) where {K,VT}
    K > 2 && (ix = MVector{K,Int}(undef))
    while length(heap) >= 2
        a = pop!(heap)
        b = pop!(heap)

        large, small = a >= b ? (a,b) : (b,a)

        subsums = ntuple(Val(K)) do i
            sums(large, i) + sums(small, K+1-i)
        end

        un = ismutable(subs(large, 1)) ? union! : union
        subsets = ntuple(Val(K)) do i
            un(subs(large, i), subs(small, K+1-i))
        end

        K<=2 && (push!(heap, VT(value(large)-value(small), subsets, subsums)); continue)

        sortperm!(ix, SVector(subsums); rev=true)
        subsums_s = subsums[ix]
        node = VT(first(subsums_s)-last(subsums_s), subsets[ix], subsums_s)
        push!(heap, node)
    end
    return pop!(heap)
end

function ldm(numbers::Union{AbstractVector,AbstractSet}, ::Val{K}, ::Type{T}) where {K,T<:Unsigned}
    ldm(numbers, Val(K), SmallBitSet{T})
end
