#
# compositions
#

export compositions, weakcompositions, compositions_cumsum, weakcompositions_cumsum

import Base: IteratorSize, HasLength, eltype, length, iterate

const CompN = 16
const CompEltype = Int8

export WeakCompositions

struct WeakCompositions
    n::Int  # may be negative (for compositions)
    k::Int  # always >= 0
end

"""
    weakcompositions_cumsum(n::Integer, k::Integer)

Return an iterator over the cumulative sums of the weak compositions of `n` of length `k`.
A weak composition of `n` of length `k` is a `k`-tuple of non-negative integers that add up to `n`.
The cumulative sum of such a composition is a vector with `k+1` elements, starting with `0` and ending with `n`.
Each vector is of type `SmallVector{$CompN,$CompEltype}`, but this may change in the future.

See also [`weakcompositions`](@ref), [`compositions_cumsum`](@ref).

# Examples
```jldoctest
julia> weakcompositions_cumsum(3, 2) |> collect
4-element Vector{SmallVector{16, Int8}}:
 [0, 0, 3]
 [0, 1, 3]
 [0, 2, 3]
 [0, 3, 3]

julia> weakcompositions_cumsum(3, 0) |> collect
SmallVector{16, Int8}[]

julia> weakcompositions_cumsum(0, 0) |> collect
1-element Vector{SmallVector{16, Int8}}:
 [0]
```
"""
function weakcompositions_cumsum(n::Integer, k::Integer)
    (n >= 0 && k >= 0) || error("arguments must be non-negative")
    k < CompN || error("only compositions into at most $(CompN-1) parts are supported")
    WeakCompositions(n, k)
end

IteratorSize(::WeakCompositions) = HasLength()

function length(c::WeakCompositions)
    if c.n+c.k-1 >= 0
        binomial(c.n+c.k-1, c.k-1)
    else #  we have c.n <= 0 && c.k == 0
        Int(iszero(c.n))
    end
end

eltype(::Type{WeakCompositions}) = SmallVector{CompN,CompEltype}

@inline function iterate(c::WeakCompositions)
    if c.n < 0 || (c.n > 0 && c.k == 0)
        nothing
    else
        v = zeros(MutableSmallVector{CompN,CompEltype}, c.k+1)
        @inbounds v[end] = c.n % CompEltype
        SmallVector(v), (v, c.k+1)
    end
end

@inline function iterate(c::WeakCompositions, (v, i)::Tuple{MutableSmallVector,Int})
    while i > 0 && @inbounds v[i] == c.n
        i -= 1
    end
    i <= 1 && return nothing
    a = @inbounds v[i] += one(CompEltype)
    for j in i+1:c.k
        @inbounds v[j] = v[i]
    end
    SmallVector(v), (v, c.k)
end

"""
    compositions_cumsum(n::Integer, k::Integer)

Return an iterator over the cumulative sums of the compositions of `n` of length `k`.
A composition of `n` of length `k` is a `k`-tuple of positive integers that add up to `n`.
The cumulative sum of such a composition is a vector with `k+1` elements, starting with `0` and ending with `n`.
Each vector is of type `SmallVector{$CompN,$CompEltype}`, but this may change in the future.

See also [`compositions`](@ref), [`weakcompositions_cumsum`](@ref).

# Examples
```jldoctest
julia> compositions_cumsum(3, 2) |> collect
2-element Vector{SmallVector{16, Int8}}:
 [0, 1, 3]
 [0, 2, 3]

julia> compositions_cumsum(3, 0) |> collect
SmallVector{16, Int8}[]

julia> compositions_cumsum(0, 0) |> collect
1-element Vector{SmallVector{16, Int8}}:
 [0]
```
"""
function compositions_cumsum(n::Integer, k::Integer)
    (n >= 0 && k >= 0) || error("arguments must be non-negative")
    k < CompN || error("only compositions into at most $(CompN-1) parts are supported")
    Generator(Fix2(composition_cumsum, k), WeakCompositions(n-k, k))
end

@inline function composition_cumsum(v::AbstractSmallVector{N,T}, k) where {N,T}
    t = ntuple(i -> ifelse(i <= k+1, T(i-1), zero(T)), Val(N))
    w = SmallVector(FixedVector(t), k+1)
    @inbounds v+w
end

eltype(::Type{<:Generator{WeakCompositions, <:Fix2{typeof(composition_cumsum)}}}) = eltype(WeakCompositions)

"""
    weakcompositions(n::Integer, k::Integer)

Return an iterator over the weak compositions of `n` of length `k`.
A weak composition of `n` of length `k` is a `k`-tuple of non-negative integers that add up to `n`.
Each composition is of type `SmallVector{$CompN,$CompEltype}`, but this may change in the future.

See also [`compositions`](@ref), [`weakcompositions_cumsum`](@ref).

# Examples
```jldoctest
julia> weakcompositions(3, 2) |> collect
4-element Vector{SmallVector{16, Int8}}:
 [0, 3]
 [1, 2]
 [2, 1]
 [3, 0]

julia> weakcompositions(3, 0) |> collect
SmallVector{16, Int8}[]

julia> weakcompositions(0, 0) |> collect
1-element Vector{SmallVector{16, Int8}}:
 0-element SmallVector{16, Int8}
```
"""
weakcompositions(n::Integer, k::Integer) = Generator(Fix2(weakcomposition, k), weakcompositions_cumsum(n, k))

@inline function weakcomposition(v::AbstractSmallVector{N,T}, k) where {N,T}
    # @inbounds first(popfirst(v)) - first(pop(v))  # too slow
    b = circshift(fixedvector(v), Val(-1)) - fixedvector(v)
    SmallVector(padtail(b, k), k)
end

eltype(::Type{<:Generator{WeakCompositions, <:Fix2{typeof(weakcomposition)}}}) = eltype(WeakCompositions)

"""
    compositions(n::Integer, k::Integer)

Return an iterator over the compositions of `n` of length `k`.
A composition of `n` of length `k` is a `k`-tuple of positive integers that add up to `n`.
Each composition is of type `SmallVector{$CompN,$CompEltype}`, but this may change in the future.

See also [`weakcompositions`](@ref), [`compositions_cumsum`](@ref).

# Examples
```jldoctest
julia> compositions(3, 2) |> collect
2-element Vector{SmallVector{16, Int8}}:
 [1, 2]
 [2, 1]

julia> compositions(3, 0) |> collect
SmallVector{16, Int8}[]

julia> compositions(0, 0)  |> collect
1-element Vector{SmallVector{16, Int8}}:
 0-element SmallVector{16, Int8}
```
"""
function compositions(n::Integer, k::Integer)
    (n >= 0 && k >= 0) || error("arguments must be non-negative")
    k < CompN || error("only compositions into at most $(CompN-1) parts are supported")
    Generator(Fix2(composition, k), WeakCompositions(n-k, k))
end

@inline function composition(v::AbstractSmallVector{N,T}, k) where {N,T}
    t = ntuple(i -> T(i <= k), Val(N))
    w = SmallVector(FixedVector(t), k)
    @inbounds weakcomposition(v, k)+w
end

eltype(::Type{<:Generator{WeakCompositions, <:Fix2{typeof(composition)}}}) = eltype(WeakCompositions)
