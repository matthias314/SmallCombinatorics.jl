#
# permutations
#

export Permutations, permutations, permutations_parity_transposition

using Base: Generator, HasEltype
import Base: length, eltype, iterate, IteratorEltype

struct Permutations
    n::Int
end

const PermN = 16
const PermEltype = Int8

"""
    permutations(n::Integer)

Return an iterator that yields all permutations of the integers from `1` to `n`.

The argument `n` must be between `0` and `$PermN`.
The identity permutation is returned first.
Each permutation is of type `SmallVector{$PermN,$PermEltype}`, but this may change in the future.

See also [`permutations_parity_transposition`](@ref).

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

const Permutations1 = Generator{Permutations, typeof(first)}

IteratorEltype(::Type{Permutations1}) = HasEltype()

eltype(::Type{Permutations1}) = SmallVector{PermN,PermEltype}

"""
    permutations_parity_transposition(n::Integer)

Return an iterator that yields all permutations `p` of the integers from `1` to `n`
together with some extra data. The first element of the tuple returned is the permutation `p`.
The second element is the parity of `p` (`false` for even and `true` for odd permutations).
The third element is a pair `(i, j)` that indicates the transposition `t` by which `p` differs
from the previously returned permutation `q`. (More precisely, the new permutations `p` is
obtained by first applying `t` and then `q`.)

The argument `n` must be between `0` and `$PermN`.
The iterator returns the identity permutation first;
in this case the transposition pair is set to `(0, 0)`.
The true transpositions `(i, j)` satisfy `i < j`.
Each permutation is of type `SmallVector{$PermN,$PermEltype}`, but this may change in the future.

See also [`permutations`](@ref).

# Examples
```jldoctest
julia> permutations_parity_transposition(3) |> collect
6-element Vector{Tuple{SmallVector{16, Int8}, Int64, Tuple{Int64, Int64}}}:
 ([1, 2, 3], 0, (0, 0))
 ([2, 1, 3], 1, (1, 2))
 ([3, 1, 2], 0, (1, 3))
 ([1, 3, 2], 1, (1, 2))
 ([2, 3, 1], 0, (1, 3))
 ([3, 2, 1], 1, (1, 2))

julia> permutations_parity_transposition(0) |> collect
1-element Vector{Tuple{SmallVector{16, Int8}, Int64, Tuple{Int64, Int64}}}:
 ([], 0, (0, 0))
```
"""
function permutations_parity_transposition(n::Integer)
    0 <= n <= PermN || error("argument must be between 0 and $PermN")
    Permutations(n)
end

length(perm::Permutations) = factorial(perm.n)

eltype(::Type{Permutations}) = Tuple{SmallVector{PermN,PermEltype},Int,NTuple{2,Int}}

# we use Heap's algorithm to iterate over all permutations

@propagate_inbounds function swap!(v::AbstractVector, i, j)
    v[i], v[j] = v[j], v[i]
    v
end

@inline function iterate(perm::Permutations)
    p = MutableSmallVector{PermN}(Base.OneTo(PermEltype(perm.n)))
    (SmallVector(p), false, (0, 0)), (p, zero(p), false)
end

@inline function iterate(perm::Permutations, (p, c, s)::Tuple{MutableSmallVector,MutableSmallVector,Bool})
    i = 1
    @inbounds while i < perm.n
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

# permutations of vectors

"""
    permutations(v::Union{AbstractFixedVector, AbstractSmallVector, PackedVector})

Return an iterator that yields all permutations of the vector `v`, whose elements are assumed to be distinct.

The vector `v` must be of length at most `$PermN`.
The elements of the iterator are of type `SmallVector{$PermN,T}` where `T` is the element type of `v`,
but this may change in future versions.
The identity permutation is returned first.

See also [`permutations(::Integer)`](@ref), [`multiset_permutations`](@ref).

# Example
```jldoctest
julia> permutations(FixedVector{3}('a':'c')) |> collect
6-element Vector{SmallVector{16, Char}}:
 ['a', 'b', 'c']
 ['b', 'a', 'c']
 ['c', 'a', 'b']
 ['a', 'c', 'b']
 ['b', 'c', 'a']
 ['c', 'b', 'a']
```
"""
permutations(v::AbstractFixedOrSmallOrPackedVector) = generator(Fix1(_inbounds_getindex, v), permutations(length(v)))

const PermutationsVector{V} = Generator{fieldtype(Permutations1, :iter), ComposedFunction{Fix1{typeof(_inbounds_getindex), V}, fieldtype(Permutations1, :f)}}

IteratorEltype(::Type{<:PermutationsVector}) = HasEltype()

eltype(::Type{PermutationsVector{V}}) where {N, T, V <: AbstractFixedOrSmallVector{N,T}} = SmallVector{N,T}
eltype(::Type{PermutationsVector{V}}) where V <: PackedVector = V

# multiset permutations

export multiset_permutations

struct MultisetPermutations{V<:SmallVector}
    v::V
end

"""
    multiset_permutations(v::AbstractSmallVector{N,T}; [sorted = false]) where {N,T}

Return an iterator over all multiset permutations of `v`, that is, all permutations
where equal elements are not distinguished. The element type `T` must have an ordering.
If `sorted` is `true`, then `v` is assumed to be sorted. The iterator yields
elements of type `SmallVector{N,T}`.

At present, the element type must satisfy `isbitstype(T)`.

See also [`permutations(::AbstractSmallVector)`](@ref), `Base.isbitstype`.

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
function multiset_permutations(v::AbstractSmallVector{N,T}; sorted = false) where {N,T}
    isbitstype(T) || error("elements must be of an isbitstype type")
    w = sorted ? SmallVector(v) : SmallVector(sort(v))
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
