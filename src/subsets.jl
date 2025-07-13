#
# subset iterators
#

export setcompositions, subsets, setcompositions_parity, setcomposition_parity

using Base: @propagate_inbounds, Generator, HasEltype
import Base: eltype, length, size, IndexStyle, getindex, iterate, IteratorEltype

struct SetCompositions{N,S}
    set::S
    ks::NTuple{N,Int}
end

"""
    setcompositions_parity(s::S, ks::Vararg{Integer,N}) where {S <: SmallBitSet, N}
    setcompositions_parity(ks::Vararg{Integer,N}) where N

In the first form, return an iterator that yields all `ks`-compositions of the set `s`
together with the parity of the permutation that puts the elements back into an increasing order.
See `setcompositions` and `setcomposition_parity` for details.
The iterator returns tuples `(t, p)`, where `t` is of type `NTuple{N, S}`
and the parity `p` is of type `Bool` where `false` means even and `true` means odd.
The partition sizes in `ks` must be non-negative and add up to `length(s)`.

In the second form the set `s` is taken to be `SmallBitSet(1:sum(ks))`.

See also [`setcompositions`](@ref), [`setcomposition_parity`](@ref).

# Examples
```jldoctest
julia> setcompositions_parity(SmallBitSet([2, 4, 5]), 1, 2) |> collect
3-element Vector{Tuple{Tuple{SmallBitSet{UInt64}, SmallBitSet{UInt64}}, Bool}}:
 ((SmallBitSet([2]), SmallBitSet([4, 5])), 0)
 ((SmallBitSet([4]), SmallBitSet([2, 5])), 1)
 ((SmallBitSet([5]), SmallBitSet([2, 4])), 0)

julia> all(s == setcomposition_parity(a, b) for ((a, b), s) in setcompositions_parity(1, 2))
true
```
"""
function setcompositions_parity(ks::Integer...)
    any(signbit, ks) && error("part sizes must be non-negative")
    sum(ks; init = 0) <= bitsize(UInt) || error("at most $(bitsize(UInt)) elements supported")
    SetCompositions(missing, ks)
end,
function setcompositions_parity(s::SmallBitSet, ks::Integer...)
    sum(ks; init = 0) == length(s) || error("part lengths must add up to size of the set")
    any(signbit, ks) && error("part sizes must be non-negative")
    SetCompositions(s, ks)
end

IteratorEltype(::Type{I}) where I <: SetCompositions = HasEltype()

eltype(::Type{SetCompositions{N,Missing}}) where N = Tuple{NTuple{N,SmallBitSet{UInt}}, Bool}
eltype(::Type{SetCompositions{N,S}}) where {N, S <: SmallBitSet} = Tuple{NTuple{N,S}, Bool}

length(sh::SetCompositions{0}) = 1

function length(sh::SetCompositions{N}) where N
    foldl(sh.ks[2:end]; init = (1, sh.ks[1])) do (p, k), l
        p*binomial(k+l, k), k+l
    end |> first
end

iterate(sh::SetCompositions{0}) = ((), false), nothing
iterate(sh::SetCompositions{0}, _) = nothing

iterate(sh::SetCompositions{1}) = ((sh.set,), false), nothing
iterate(sh::SetCompositions{1,Missing}) = ((SmallBitSet(1:sh.ks[1]),), false), nothing
iterate(sh::SetCompositions{1}, _) = nothing

@inline iterate(sh::SetCompositions{2}) = any(signbit, sh.ks) ? nothing : _iterate(sh)

@inline function _iterate(sh::SetCompositions{2,S}; signint = UInt(0)) where S
    k, l = sh.ks
    U = S == Missing ? UInt : typeof(bits(sh.set))
    mask = U(1) << k - U(1)
    lastmask = mask << l

    set = S == Missing ? @inbounds(SmallBitSet{U}(1:k+l)) : sh.set
    part1 = _SmallBitSet(S == Missing ? mask : pdep(mask, bits(set)))
    part2 = symdiff(part1, set)
    signbit = isodd(signint)
    state = (; mask, lastmask, signint, set, part2)
    ((part1, part2), signbit), (state,)
end

@inline function iterate(sh::SetCompositions{2,S}, (state,)) where S
    (; mask, lastmask, signint, set) = state

    # see also https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
    # and https://discourse.julialang.org/t/faster-way-to-find-all-bit-arrays-of-weight-n/113658/12
    mask == lastmask && return nothing
    p = mask + blsi(mask)
    t = trailing_zeros(mask)
    q = unsafe_lshr(mask ⊻ p, t) >>> 2
    # q = unsafe_lshr(blsmsk(p), t) >>> 2
    # t+2 can be the bit size of mask, so we can't use unsafe_lshr with t+2
    mask = p | q
    signint ⊻= ~(t & count_ones(q))

    part1 = _SmallBitSet(S == Missing ? mask : pdep(mask, bits(set)))
    part2 = symdiff(part1, set)
    signbit = isodd(signint)
    state = (; mask, lastmask, signint, set, part2)
    ((part1, part2), signbit), (state,)
end

@inline iterate(sh::SetCompositions) = _iterate(sh)

@inline function _iterate(sh::SetCompositions{N}) where N
    sh2 = SetCompositions(sh.set, (sh.ks[1]+sh.ks[2], sh.ks[3:end]...))
    ((set1, _), _), states_rest = _iterate(sh2)
    ((part1, part2), _), (state1,) = _iterate(SetCompositions(set1, sh.ks[1:2]))
    states = (state1, states_rest...)
    parts = (part1, map(state -> state.part2, states)...)
    signbit = false
    (parts, signbit), states
end

@inline function iterate(sh::SetCompositions{N}, states) where N
    ts1 = iterate(SetCompositions(states[1].set, sh.ks[1:2]), (states[1],))
    if ts1 === nothing
        sh_rest = SetCompositions(states[2].set, (sh.ks[1]+sh.ks[2], sh.ks[3:end]...))
        ts_rest = iterate(sh_rest, states[2:end])
        ts_rest === nothing && return nothing
        ((set1, _), _), states_rest = ts_rest
        ((part1, _), signbit), (state1,) = _iterate(SetCompositions(set1, sh.ks[1:2]); signint = states_rest[1].signint)
        states = (state1, states_rest...)
    else
        ((part1, _), signbit), (state1,) = ts1
        states = (state1, states[2:end]...)
    end
    parts = (part1, map(state -> state.part2, states)...)
    (parts, signbit), states
end

"""
    setcomposition_parity(ss::SmallBitSet...) -> Bool

Return `true` if an odd number of transpositions is needed to transform the elements of the
sets `ss` into an increasing sequence, and `false` otherwise. The sets are considered as
increasing sequences and assumed to be disjoint.

See also [`setcompositions_parity`](@ref).

# Examples
```
julia> s, t, u = SmallBitSet([2, 3, 8]), SmallBitSet([1, 4, 6]), SmallBitSet([5, 7]);

julia> setcomposition_parity(s, t), setcomposition_parity(s, t, u)
(true, false)
```
"""
setcomposition_parity(ss::Vararg{SmallBitSet,N}) where N =
    setcomposition_parity(ss[N-1], ss[N]) ⊻ (@inline setcomposition_parity(ss[1:N-2]..., ss[N-1] ∪ ss[N]))

setcomposition_parity() = false
setcomposition_parity(::SmallBitSet) = false

function setcomposition_parity(s::SmallBitSet, t::SmallBitSet)
    m = bits(s)
    p = zero(m)
    while !iszero(m)
        # p ⊻= blsi(m)-one(m)
        p ⊻= blsmsk(m)
        m = blsr(m)
    end
    isodd(count_ones(p & bits(t)))
end

"""
    setcompositions(s::S, ks::Vararg{Integer,N}) where {S <: SmallBitSet, N}
    setcompositions(ks::Vararg{Integer,N}) where N

In the first form, return an iterator that yields all `ks`-compositions of the set `s`, that is,
all ordered partitions of `s` into `N` sets of size `ks[1]` to `ks[N]`, respectively. The element type
is `NTuple{N, S}`. The partition sizes in `ks` must be non-negative and add up to `length(s)`.

In the second form the set `s` is taken to be `SmallBitSet(1:sum(ks))`.
This gives an iterator over all set compositions of the integer `sum(ks)`.

See also [`subsets`](@ref subsets(::SmallBitSet, ::Integer)),
[`setcompositions_parity`](@ref setcompositions_parity(::Vararg{Integer,N}) where N).

# Examples
```jldoctest
julia> setcompositions(SmallBitSet([2, 4, 5]), 1, 2) |> collect
3-element Vector{Tuple{SmallBitSet{UInt64}, SmallBitSet{UInt64}}}:
 (SmallBitSet([2]), SmallBitSet([4, 5]))
 (SmallBitSet([4]), SmallBitSet([2, 5]))
 (SmallBitSet([5]), SmallBitSet([2, 4]))

julia> setcompositions(1, 1, 1) |> collect
6-element Vector{Tuple{SmallBitSet{UInt64}, SmallBitSet{UInt64}, SmallBitSet{UInt64}}}:
 (SmallBitSet([1]), SmallBitSet([2]), SmallBitSet([3]))
 (SmallBitSet([2]), SmallBitSet([1]), SmallBitSet([3]))
 (SmallBitSet([1]), SmallBitSet([3]), SmallBitSet([2]))
 (SmallBitSet([3]), SmallBitSet([1]), SmallBitSet([2]))
 (SmallBitSet([2]), SmallBitSet([3]), SmallBitSet([1]))
 (SmallBitSet([3]), SmallBitSet([2]), SmallBitSet([1]))

julia> setcompositions(SmallBitSet([2, 4, 5]), 1, 0, 2) |> collect
3-element Vector{Tuple{SmallBitSet{UInt64}, SmallBitSet{UInt64}, SmallBitSet{UInt64}}}:
 (SmallBitSet([2]), SmallBitSet([]), SmallBitSet([4, 5]))
 (SmallBitSet([4]), SmallBitSet([]), SmallBitSet([2, 5]))
 (SmallBitSet([5]), SmallBitSet([]), SmallBitSet([2, 4]))

julia> setcompositions(SmallBitSet()) |> collect
1-element Vector{Tuple{}}:
 ()
```
"""
setcompositions(args...) = Generator(first, setcompositions_parity(args...))

IteratorEltype(::Type{Generator{I, typeof(first)}}) where I <: SetCompositions = HasEltype()

eltype(g::Type{Generator{I, typeof(first)}}) where I <: SetCompositions = fieldtype(eltype(I), 1)

struct Subsets{T,S<:SmallBitSet} <: AbstractVector{S}
    set::T
    length::Int
end

"""
    subsets(s::S) where S <: SmallBitSet -> AbstractVector{S}
    subsets(n::Integer) -> AbstractVector{SmallBitSet{UInt}}

In the first form, return a vector of length `2^length(s)` whose elements are the subsets of the set `s`.

In the second form the set `s` is taken to be `SmallBitSet(1:n)`.

See also [`subsets(::Integer, ::Integer)`](@ref).

# Examples
```jldoctest
julia> subsets(SmallBitSet{UInt8}([3, 5])) |> collect
4-element Vector{SmallBitSet{UInt8}}:
 SmallBitSet([])
 SmallBitSet([3])
 SmallBitSet([5])
 SmallBitSet([3, 5])

julia> subsets(2) |> collect
4-element Vector{SmallBitSet{UInt64}}:
 SmallBitSet([])
 SmallBitSet([1])
 SmallBitSet([2])
 SmallBitSet([1, 2])

julia> subsets(2)[2]
SmallBitSet{UInt64} with 1 element:
  1
```
"""
function subsets(n::T) where T <: Integer
    n >= 0 || error("argument must be non-negative")
    n <= bitsize(UInt)-2 || error("at most $(bitsize(UInt)-2) elements supported")
    Subsets{T,SmallBitSet{UInt}}(n, n >= 0 ? unsafe_shl(1, n) : 0)
end,
function subsets(s::S) where {U <: Unsigned, S <: SmallBitSet{U}}
    bitsize(U) < bitsize(UInt) || length(s) <= bitsize(UInt)-2 ||
        error("at most $(bitsize(UInt)-2) elements supported")
    Subsets{S,S}(s, unsafe_shl(1, length(s)))
end

show(io::IO, ss::Subsets) = print(io, "Subsets(", ss.set, ')')
show(io::IO, ::MIME"text/plain", ss::Subsets) = print(io, "Subsets(", ss.set, ')')

size(ss::Subsets) = (ss.length,)

IndexStyle(::Type{<:Subsets}) = IndexLinear()

@inline function getindex(ss::Subsets{<:Integer}, i::Int)
    @boundscheck checkbounds(ss, i)
    _SmallBitSet((i-1) % UInt)
end

@inline function getindex(ss::Subsets{<:SmallBitSet}, i::Int)
    @boundscheck checkbounds(ss, i)
    _SmallBitSet(pdep((i-1) % UInt, bits(ss.set)))
end

"""
    subsets(s::Union{SmallBitSet, AbstractSmallSet}, k::Integer)
    subsets(n::Integer, k::Integer)

In the first form, return an iterator that yields all `k`-element subsets of the set `s`.
The element type of the iterator is a `SmallBitSet` or `SmallSet`.
If `k` is negative or larger than `length(s)`, then the iterator is empty.

In the second form the set `s` is taken to be `SmallBitSet(1:n)`.

See also [`subsets(::Integer)`](@ref),
[`combinations`](@ref combinations(::Integer, ::Integer)),
[`setcompositions_parity`](@ref setcompositions_parity(::Vararg{Integer,N}) where N).

# Example
```jldoctest
julia> subsets(SmallBitSet{UInt8}(2:2:8), 3) |> collect
4-element Vector{SmallBitSet{UInt8}}:
 SmallBitSet([2, 4, 6])
 SmallBitSet([2, 4, 8])
 SmallBitSet([2, 6, 8])
 SmallBitSet([4, 6, 8])

julia> subsets(3, 2) |> collect
3-element Vector{SmallBitSet{UInt64}}:
 SmallBitSet([1, 2])
 SmallBitSet([1, 3])
 SmallBitSet([2, 3])

julia> subsets(MutableSmallSet{4}('a':'c'), 2) |> collect
3-element Vector{SmallSet{4, Char}}:
 SmallSet{4}(['a', 'b'])
 SmallSet{4}(['a', 'c'])
 SmallSet{4}(['b', 'c'])

julia> subsets(3, 4) |> collect
SmallBitSet{UInt64}[]
```
"""
subsets(::Union{SmallBitSet, AbstractSmallSet, Integer}, ::Integer)

function subsets(n::Integer, k::Integer)
    n >= 0 || error("first argument must be non-negative")
    n <= bitsize(UInt) || error("at most $(bitsize(UInt)) elements supported")
    Generator(first∘first, SetCompositions(missing, (k, n-k)))
end

function subsets(s::SmallBitSet, k::Integer)
    Generator(first∘first, SetCompositions(s, (k, length(s)-k)))
end

const Subsets2{S} = Generator{SetCompositions{2,S}, typeof(first∘first)}
const Subsets2Int = Subsets2{Missing}

IteratorEltype(::Type{<:Subsets2}) = HasEltype()

eltype(::Type{Subsets2Int}) = SmallBitSet{UInt}
eltype(::Type{Subsets2{S}}) where S <: SmallBitSet = S

subsets(s::AbstractSmallSet, k::Integer) = combinations(s, k)

# combinations

export combinations

"""
    combinations(n::Integer, k::Integer)
    combinations(s::Union{SmallBitSet, AbstractSmallSet}, k::Integer)
    combinations(s::Union{AbstractFixedVector, AbstractSmallVector}, k::Integer)
    combinations(s::PackedVector, k::Integer)

Return an iterator that yields all combinations of `k` elements from the given collection,
whose elements are assumed to be distinct.
If the first argument is an integer `n`, then the collection is taken to be `SmallBitSet(1:n)`.
If `k` is negative or exceeds the length of the collection, then the iterator is empty.

With an integer or a set as first argument, `combinations` is the same as `subsets`.
If the collection is an `AbstractFixedVector{N,T}` or `AbstractSmallVector{N,T}`, then the
element type of the iterator is `SmallVector{N,T}`. If the collection is a `PackedVector`,
then the elements of the iterator are of the same type.

The version with an integer or a `SmallBitSet` as first argument is the fastest.

See also [`subsets`](@ref subsets(::Integer, ::Integer)).

# Example
```jldoctest
julia> combinations(MutableFixedVector{4}('a':'d'), 2) |> collect
6-element Vector{SmallVector{4, Char}}:
 ['a', 'b']
 ['a', 'c']
 ['b', 'c']
 ['a', 'd']
 ['b', 'd']
 ['c', 'd']
```
"""
combinations(::Union{Integer, SmallBitSet, AbstractSmallSet, AbstractFixedOrSmallOrPackedVector}, ::Integer)

combinations(n::Integer, k::Integer) = subsets(n, k)
combinations(s::SmallBitSet, k::Integer) = subsets(s, k)

_inbounds_getindex(v::AbstractFixedOrSmallOrPackedVector, ii) = @inbounds v[ii]
_inbounds_getindex(s::AbstractSmallSet{N,T}, ii) where {N,T} = SmallSet{N,T}(@inbounds values(s)[ii]; unique = true)

combinations(c::Union{AbstractFixedOrSmallOrPackedVector, AbstractSmallSet}, k::Integer) = generator(Fix1(_inbounds_getindex, c), subsets(length(c), k))

const Combinations2{C} = Generator{fieldtype(Subsets2Int, :iter), ComposedFunction{Fix1{typeof(_inbounds_getindex), C}, fieldtype(Subsets2Int, :f)}}

IteratorEltype(::Type{<:Combinations2}) = HasEltype()

eltype(::Type{Combinations2{V}}) where {N, T, V <: AbstractFixedOrSmallVector{N,T}} = SmallVector{N,T}
eltype(::Type{Combinations2{V}}) where V <: PackedVector = V
eltype(::Type{Combinations2{S}}) where {N, T, S <: AbstractSmallSet{N,T}} = SmallSet{N,T}
