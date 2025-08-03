#
# partitions
#

export partitions

import Base: IteratorSize, eltype, iterate

const PartN = 64

const PartEltype = Int8

struct Partitions{V<:AbstractSmallVector}
    n::Int
end

"""
    partitions(n::Integer)
    partitions(::Type{V}, n::Integer) where {N, T <: Base.BitInteger, V <: AbstractSmallVector{N,T}}

Return an iterator over the partitions of `n`.
A partition of `n` is a weakly decreasing sequence of positive integers that add up to `n`.
Each partition is of type `V <: AbstractSmallVector` if this type is specified.
Otherwise, `V` is taken to be `SmallVector{$PartN,$PartEltype}`; this may change in the future.

See also [`partitions(::Integer, ::Integer)`](@ref), `Base.BitInteger`.

# Examples
```jldoctest
julia> partitions(3) |> collect
3-element Vector{SmallVector{64, Int8}}:
 [3]
 [2, 1]
 [1, 1, 1]

julia> partitions(SmallVector{8,UInt8}, 3) |> collect
3-element Vector{SmallVector{8, UInt8}}:
 [0x03]
 [0x02, 0x01]
 [0x01, 0x01, 0x01]

julia> partitions(0) |> collect
1-element Vector{SmallVector{64, Int8}}:
 0-element SmallVector{64, Int8}
```
"""
partitions(::Integer),
partitions(::Type{<:AbstractSmallVector{<:Any,<:Base.BitInteger}}, ::Integer)

partitions(n::Integer) = partitions(SmallVector{PartN,PartEltype}, n)

function partitions(::Type{V}, n::Integer) where {N, T <: Base.BitInteger, V <: AbstractSmallVector{N,T}}
    M = min(N, typemax(T))
    0 <= n <= M || error("argument must be between 0 and $M for $V")
    Partitions{V}(n)
end

IteratorSize(::Type{<:Partitions}) = Base.SizeUnknown()

eltype(::Type{Partitions{V}}) where V = V

@inline function iterate(p::Partitions{V}) where {N, T, V <: AbstractSmallVector{N,T}}
    v = MutableSmallVector{N,T}()
    # creating two separate MutableSmallVector below would lead to allocations
    if iszero(p.n)
        V(v), (v, 0, zero(T))
    else
        @inbounds push!(v, T(p.n))
        V(v), (v, 1-isone(p.n), T(isone(p.n)))
    end
end

@inline function iterate(p::Partitions{V}, (v, i, s)::Tuple{MutableSmallVector,Int,T}) where {N, T, V <: AbstractSmallVector{N,T}}
    @inbounds if i == 0
        nothing
    elseif v[i] == 2
        v[i] = 1
        push!(v, one(T))
        V(v), (v, i-1, s+T(2))
    else
        c = (v[i] -= one(T))
        s += one(T)
        while s >= c
            i += 1
            v[i] = c
            s -= c
        end
        if s == 0
            resize!(v, i)
            V(v), (v, length(v), zero(T))
        else
            resize!(v, i+1)
            v[i+1] = s
            V(v), (v, length(v)-isone(s), T(isone(s)))
        end
    end
end

struct Partitions2{V<:AbstractSmallVector}
    n::Int
    k::Int
end

"""
    partitions(n::Integer, k::Integer)
    partitions(::Type{V}, n::Integer, k::Integer) where {N, T <: Base.BitInteger, V <: AbstractSmallVector{N,T}}

Return an iterator over the partitions of `n` into `k` parts.
A partition of `n` is a weakly decreasing sequence of positive integers that add up to `n`.
Each partition is of type `V <: AbstractSmallVector` if this type is specified.
Otherwise, `V` is taken to be `SmallVector{$PartN,$PartEltype}`; this may change in the future.

See also [`partitions(::Integer)`](@ref), `Base.BitInteger`.

# Examples
```jldoctest
julia> partitions(7, 3) |> collect
4-element Vector{SmallVector{64, Int8}}:
 [5, 1, 1]
 [4, 2, 1]
 [3, 3, 1]
 [3, 2, 2]

julia> partitions(SmallVector{8,UInt8}, 7, 3) |> collect
4-element Vector{SmallVector{8, UInt8}}:
 [0x05, 0x01, 0x01]
 [0x04, 0x02, 0x01]
 [0x03, 0x03, 0x01]
 [0x03, 0x02, 0x02]

julia> partitions(7, 0) |> collect
SmallVector{64, Int8}[]

julia> partitions(0, 0) |> collect
1-element Vector{SmallVector{64, Int8}}:
 0-element SmallVector{64, Int8}
```
"""
partitions(::Integer, ::Integer),
partitions(::Type{<:AbstractSmallVector{<:Any,<:Base.BitInteger}}, ::Integer, ::Integer)

partitions(n::Integer, k::Integer) = partitions(SmallVector{PartN,PartEltype}, n, k)

function partitions(::Type{V}, n::Integer, k::Integer) where {N, T <: Base.BitInteger, V <: AbstractSmallVector{N,T}}
    (n >= 0 && k >= 0) || error("integer arguments must be non-negative")
    n <= typemax(T) || error("only partitions of integers up to $(typemax(T)) are supported for $V")
    k <= N || error("partitions into at most $N parts are supported for $V")
    Partitions2{V}(n, k)
end

IteratorSize(::Type{<:Partitions2}) = Base.SizeUnknown()

eltype(::Type{Partitions2{V}}) where V = V

@inline function iterate(p::Partitions2{V}) where {N, T, V <: AbstractSmallVector{N,T}}
    (p.n < p.k || p.n > p.k == 0) && return nothing
    v = MutableSmallVector{N,T}(undef, p.k)
    for i in eachindex(v)
        @inbounds v[i] = (i == 1 ? p.n-p.k+1 : 1) % T
    end
    V(v), v
end

@inline function iterate(p::Partitions2{V}, v) where {N, T, V <: AbstractSmallVector{N,T}}
    @inbounds begin
        local c
        if p.k == 0 || (c = v[1] - one(T); c <= v[p.k])
            return nothing
        end
        i = findfirst(<(c), v)::Int
        c = (v[i] += one(T))
        # unsafe_copyto!(v, map(Fix2(min, c), v; style = RigidStyle()))  # allocates
        unsafe_copyto!(v, min.(fixedvector(v), c))
        v[1] += p.n % T - sum_fast(v)
    end
    V(v), v
end
