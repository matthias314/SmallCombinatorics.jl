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
