#
# permutations
#

export permutations, permutations_parity_transposition

using Base: Generator, HasEltype
import Base: length, eltype, iterate, IteratorEltype

struct PermutationsParityTransposition{V<:SmallVector}
    v::V
end

const PermN = 16
const PermEltype = Int8

"""
    permutations(n::Integer)

Return an iterator that yields all permutations of the integers from `1` to `n`.
The argument `n` must be between `0` and `$PermN`.

This is a short form for `permutations(v)` with `v = SmallVector{$PermN,$PermEltype}(1:n))`.
(Capacity and element type of `v` may change in the future.)

See also [`permutations(::SmallVector)`](@ref), [`permutations_parity_transposition`](@ref).

# Examples
```jldoctest
julia> permutations(3) |> collect
6-element Vector{SmallVector{16, Int8}}:
 [1, 2, 3]
 [2, 1, 3]
 [3, 1, 2]
 [1, 3, 2]
 [2, 3, 1]
 [3, 2, 1]

julia> permutations(0) |> collect
1-element Vector{SmallVector{16, Int8}}:
 0-element SmallVector{16, Int8}
```
"""
permutations(n::Integer) = Generator(first, permutations_parity_transposition(n))

const Permutations{V} = Generator{PermutationsParityTransposition{V}, typeof(first)}

IteratorEltype(::Type{Permutations}) = HasEltype()

eltype(::Type{Permutations{V}}) where V = V

"""
    permutations(v::Union{AbstractSmallVector, AbstractFixedVector, PackedVector})

Return an iterator that yields all permutations of the vector `v`, whose elements are assumed to be distinct.

At present, the elements of the iterator are of type `SmallVector{N,T}` where `N` is the capacity of `v` and `T` the element type.
The identity permutation is returned first.

See also [`permutations(::Integer)`](@ref), [`permutations_parity_transposition`](@ref), [`multiset_permutations`](@ref).

# Example
```jldoctest
julia> permutations(FixedVector{3}('a':'c')) |> collect
6-element Vector{SmallVector{3, Char}}:
 ['a', 'b', 'c']
 ['b', 'a', 'c']
 ['c', 'a', 'b']
 ['a', 'c', 'b']
 ['b', 'c', 'a']
 ['c', 'b', 'a']

julia> permutations(SmallVector{4,Symbol}()) |> collect
1-element Vector{SmallVector{4, Symbol}}:
 0-element SmallVector{4, Symbol}
```
"""
function permutations(v::AbstractFixedOrSmallOrPackedVector)
    capacity(v) <= typemax(PermEltype) || length(v) <= typemax(PermEltype) || error("vector length can be at most ", typemax(PermEltype))
    if isbitstype(eltype(v))
        Generator(first, PermutationsParityTransposition(SmallVector(v)))
    else
        N = capacity(v)
        w = SmallVector{N}(Base.OneTo(PermEltype(length(v))))
        generator(Fix1(_inbounds_getindex, v), permutations(w))
    end
end

const PermutationsGenerator{V,W} =
    Generator{fieldtype(Permutations{W}, :iter), ComposedFunction{Fix1{typeof(_inbounds_getindex), V}, fieldtype(Permutations{W}, :f)}}

IteratorEltype(::Type{<:PermutationsGenerator}) = HasEltype()

eltype(::Type{PermutationsGenerator{V,SmallVector{N,PermEltype}}}) where {V,N} = SmallVector{N,eltype(V)}

"""
    permutations_parity_transposition(n::Integer)

This is a short form for `permutations_parity_transposition(v)` with `v = SmallVector{$PermN,$PermEltype}(1:n))`.
(Capacity and element type of `v` may change in the future.)
The argument `n` must be between `0` and `$PermN`.

See also [`permutations`](@ref), [`permutations_parity_transposition(n::SmallVector)`](@ref).

# Examples
```jldoctest
julia> permutations_parity_transposition(3) |> collect
6-element Vector{Tuple{SmallVector{16, Int8}, Bool, Tuple{Int64, Int64}}}:
 ([1, 2, 3], 0, (0, 0))
 ([2, 1, 3], 1, (1, 2))
 ([3, 1, 2], 0, (1, 3))
 ([1, 3, 2], 1, (1, 2))
 ([2, 3, 1], 0, (1, 3))
 ([3, 2, 1], 1, (1, 2))

julia> permutations_parity_transposition(0) |> collect
1-element Vector{Tuple{SmallVector{16, Int8}, Bool, Tuple{Int64, Int64}}}:
 ([], 0, (0, 0))
```
"""
function permutations_parity_transposition(n::Integer)
    0 <= n <= PermN || error("argument must be between 0 and $PermN")
    permutations_parity_transposition(SmallVector{PermN}(Base.OneTo(PermEltype(n))))
end

"""
    permutations_parity_transposition(v::Union{AbstractSmallVector, AbstractFixedVector, PackedVector)

Return an iterator that yields all permutations `p` of the elements of `v`
together with some extra data. The first element of the tuple returned is the permutation `p`.
The second element is the parity of `p` (`false` for even and `true` for odd permutations).
The third element is a pair `(i, j)` that indicates the transposition `t` by which `p` differs
from the previously returned permutation `q`. (More precisely, the new permutations `p` is
obtained by first applying `t` and then `q`.)

The iterator returns the identity permutation first;
in this case the transposition pair is set to `(0, 0)`.
The true transpositions `(i, j)` satisfy `i < j`.
At present, each permutation is of type `SmallVector{N,T}` where `N` is the capacity of `v` and `T` the element type.

See also [`permutations`](@ref), [`permutations_parity_transposition(n::Integer)`](@ref).

Examples
```jldoctest
julia> v = SmallVector{4}('a':'c');

julia> permutations_parity_transposition(v) |> collect
6-element Vector{Tuple{SmallVector{4, Char}, Bool, Tuple{Int64, Int64}}}:
 (['a', 'b', 'c'], 0, (0, 0))
 (['b', 'a', 'c'], 1, (1, 2))
 (['c', 'a', 'b'], 0, (1, 3))
 (['a', 'c', 'b'], 1, (1, 2))
 (['b', 'c', 'a'], 0, (1, 3))
 (['c', 'b', 'a'], 1, (1, 2))

julia> v = SmallVector{4,String}();

julia> permutations_parity_transposition(v) |> collect
1-element Vector{Tuple{SmallVector{4, String}, Bool, Tuple{Int64, Int64}}}:
 ([], 0, (0, 0))
```
"""
function permutations_parity_transposition(v::AbstractFixedOrSmallOrPackedVector)
    capacity(v) <= typemax(PermEltype) || length(v) <= typemax(PermEltype) || error("vector length can be at most ", typemax(PermEltype))
    if isbitstype(eltype(v))
        PermutationsParityTransposition(SmallVector(v))
    else
        N = capacity(v)
        w = SmallVector{N}(Base.OneTo(PermEltype(length(v))))
        Generator(Fix1(_inbounds_getindex_first, v), PermutationsParityTransposition(w))
    end
end

_inbounds_getindex_first(v::AbstractFixedOrSmallVector, t::Tuple) = (_inbounds_getindex(v, t[1]), t[2:end]...)

length(perm::PermutationsParityTransposition) = factorial(length(perm.v))

eltype(::Type{PermutationsParityTransposition{V}}) where V = Tuple{V, Bool, NTuple{2, Int}}

const PermutationsParityTranspositionGenerator{V,W} =
    Generator{PermutationsParityTransposition{W}, Fix1{typeof(_inbounds_getindex_first), V}}

IteratorEltype(::Type{<:PermutationsParityTranspositionGenerator}) = HasEltype()

eltype(::Type{PermutationsParityTranspositionGenerator{V,SmallVector{N,PermEltype}}}) where {V,N} = Tuple{SmallVector{N,eltype(V)}, Bool, NTuple{2, Int}}

# we use Heap's algorithm to iterate over all permutations

@propagate_inbounds function swap!(v::AbstractVector, i, j)
    v[i], v[j] = v[j], v[i]
    v
end

@inline function iterate(perm::PermutationsParityTransposition{SmallVector{N,T}}) where {N,T}
    (perm.v, false, (0, 0)), (MutableSmallVector(perm.v), zeros(MutableSmallVector{N,PermEltype}, length(perm.v)), false)
end

@inline function iterate(perm::PermutationsParityTransposition, (p, c, s)::Tuple{MutableSmallVector,MutableSmallVector,Bool})
    i = 1
    @inbounds while i < length(perm.v)
        if c[i] < i
            t = iseven(i) ? (swap!(p, 1, i+1); (1, i+1)) : (swap!(p, c[i]+1, i+1); (c[i]+1, i+1))
            c[i] += one(PermEltype)
            return (SmallVector(p), !s, t), (p, c, !s)
        else
            c[i] = 0
            i += 1
        end
    end
    nothing
end

# multiset permutations

export multiset_permutations

struct MultisetPermutations{V<:SmallVector}
    v::V
end

"""
    multiset_permutations(v::Union{AbstractSmallVector, AbstractFixedVector, PackedVector}; [sorted = false])

Return an iterator over all multiset permutations of `v`, that is, all permutations
where equal elements are not distinguished. The element type `T` must have an ordering.
If `sorted` is `true`, then `v` is assumed to be sorted.

At present, the element type must satisfy `isbitstype(T)`,
and the iterator yields elements of type `SmallVector{N,T}` where `N` is the capacity of `v`.

See also [`permutations`](@ref), `Base.isbitstype`.

# Example
```jldoctest
julia> v = SmallVector{8,Int8}([1, 2, 1]);

julia> multiset_permutations(v) |> collect
3-element Vector{SmallVector{8, Int8}}:
 [1, 1, 2]
 [1, 2, 1]
 [2, 1, 1]
```
"""
function multiset_permutations(v::AbstractFixedOrSmallOrPackedVector; sorted = false)
    isbitstype(eltype(v)) || error("elements must be of an isbitstype type")
    w = sorted ? SmallVector(v) : SmallVector(sort(SmallVector(v)))
    MultisetPermutations(w)
end

@inline function iterate(mp::MultisetPermutations)
    mp.v, MutableSmallVector(mp.v)
end

@inline function iterate(mp::MultisetPermutations, w::MutableSmallVector)
    # Pandita's algorithm
    isempty(w) && return nothing
    @inbounds u, _ = pop(w)
    @inbounds v, _ = popfirst(w)
    k = findlast(map(>, v, u))
    k === nothing && return nothing
    l = findlast(>(@inbounds w[k]), w)::Int
    @inbounds swap!(w, k, l)
    @inbounds reverse!(w, k+1)
    SmallVector(w), w
end

IteratorEltype(::Type{<:MultisetPermutations}) = HasEltype()

eltype(mp::MultisetPermutations{V}) where V = V

function mp_length(::Type{T}, v::AbstractVector) where T
    n = i = T(1)
    for k in T(2):(length(v) % T)
        @inbounds i = v[k] == v[k-1] ? i+T(1) : T(1)
        n = div(n*k, i)
    end
    n
end

function length(mp::MultisetPermutations)
    # 32-bit div is faster, and factorial(12) <= typemax(UInt32)
    length(mp.v) <= 12 ? mp_length(UInt32, mp.v) % Int : mp_length(Int, mp.v)
end
