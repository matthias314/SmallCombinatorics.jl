@testset "subsets(n,k)" begin
    for n in [-1, 0, 1, 2, 10], k in [-1, 0, 1, n-1, n, n+1]
        if n < 0
            @test_throws Exception subsets(n, k)
            continue
        end
        ss = @inferred subsets(n, k)
        if 0 <= k <= n
            @test_inferred length(ss) binomial(n, k)
        else
            @test_inferred length(ss) 0
        end
        @test eltype(ss) == SmallBitSet{UInt}
        ssv = @inferred collect(ss)
        @test length(ssv) == length(ss) == length(unique(ssv))
        length(ss) == 0 && continue
        @test unique(map(length, ssv)) == [k]
    end

    for U in unsigned_types, a in map(SmallBitSet{U}, [Int[], [3], [bitsize(U)-3, bitsize(U)], [2, 4, 6, 7]])
        n = length(a)
        for k in [-1, 0, 1, n-1, n, n+1]
            ss = @inferred subsets(a, k)
            if 0 <= k <= n
                @test_inferred length(ss) binomial(n, k)
            else
                @test_inferred length(ss) 0
            end
            @test eltype(ss) == SmallBitSet{U}
            ssv = @inferred collect(ss)
            @test length(ssv) == length(ss) == length(unique(ssv))
            length(ss) == 0 && continue
            @test unique(map(length, ssv)) == [k]
        end
    end

    @test_inferred eltype(typeof(subsets(5, 2))) SmallBitSet{UInt}
    @test_inferred eltype(typeof(subsets(SmallBitSet{UInt8}(1:5), 2))) SmallBitSet{UInt8}
end

@testset "subsets(n)" begin
    for n in [-1, 0, 1, 2, 10]
        if n < 0
            @test_throws Exception subsets(n)
            continue
        end
        ss = subsets(n)
        @test_inferred length(ss) 2^n
        @test eltype(ss) == SmallBitSet{UInt}
        ssv = @inferred collect(ss)
        @test length(ssv) == length(ss) == length(unique(ssv))
        @test_throws BoundsError ss[firstindex(ss)-1]
        @test_throws BoundsError ss[lastindex(ss)+1]
    end

    for U in unsigned_types, a in map(SmallBitSet{U}, [Int[], [3], [bitsize(U)-3, bitsize(U)], [2, 4, 6, 7]])
        ss = subsets(a)
        @test_inferred length(ss) 2^length(a)
        @test eltype(ss) == SmallBitSet{U}
        ssv = @inferred collect(ss)
        @test length(ssv) == length(ss) == length(unique(ssv))
        @test_throws BoundsError ss[firstindex(ss)-1]
        @test_throws BoundsError ss[lastindex(ss)+1]
    end
end

@testset "setcompositions_parity" begin
    function test_setcompositions_parity(ks)
        N = length(ks)
        sh = @inferred setcompositions_parity(ks...)
        a = SmallBitSet(1:sum(ks; init = 0))
        @test @inferred(length(sh)) == factorial(big(sum(ks; init = 0))) ÷ prod(map(factorial∘big, ks); init = 1)
        @test @inferred(eltype(sh)) == Tuple{NTuple{N, SmallBitSet{UInt}}, Bool}
        @test all(map(length, t) == ks && s isa Bool &&
            (isempty(t) ? isempty(a) : (union(t...) == a)) &&
            setcomposition_parity(t...) == s for (t, s) in sh)
        @test allunique(sh)
    end

    function test_setcompositions_parity(a::S, ks::NTuple{N,Int}) where {S <: SmallBitSet, N}
        test_setcompositions_parity(ks)
        sh = @inferred setcompositions_parity(a, ks...)
        @test @inferred(length(sh)) == factorial(big(sum(ks; init = 0))) ÷ prod(map(factorial∘big, ks); init = 1)
        @test @inferred(eltype(sh)) == Tuple{NTuple{N, S}, Bool}
        @test all(t isa NTuple{N, S} && s isa Bool &&
            map(length, t) == ks &&
            (isempty(t) ? isempty(a) : (union(t...) == a)) &&
            setcomposition_parity(t...) == s for (t, s) in sh)
        @test allunique(sh)
    end

    for U in unsigned_types, (v, ks) in [(Int[], ()), (Int[], (0,)), (Int[], (0, 0)),
                (bitsize(U)-4:2:bitsize(U), (1, 1, 1)),
                (3:2:11, (5,)), (3:2:11, (2, 3)),  (3:2:11, (0, 2, 3)),  (3:2:11, (2, 0, 3)),  (3:2:11, (2, 3, 0)),
                (20:2:38, (2, 3, 2, 3)), (20:2:38, (1, 4, 0, 2, 3))]
        maximum(v; init = 0) <= bitsize(U) || continue
        a = SmallBitSet{U}(v)
        test_setcompositions_parity(a, ks)
    end

    @test_throws Exception setcompositions_parity(-1, 2)
    @test_throws Exception setcompositions_parity(bitsize(UInt)-1, 2)
    for U in unsigned_types
        @test_throws Exception setcompositions_parity(SmallBitSet{U}(2:2:6))
        @test_throws Exception setcompositions_parity(SmallBitSet{U}(2:2:6), -1, 2, 2)
        @test_throws Exception setcompositions_parity(SmallBitSet{U}(2:2:6), 3, 4)
        @test (setcompositions_parity(SmallBitSet{U}(1:bitsize(U)), bitsize(U)-2, 2); true)
    end

    # check that unsafe_lshr in iterate for SetCompositions is safe
    @test collect(subsets(bitsize(UInt), 1)) == [SmallBitSet((k,)) for k in 1:bitsize(UInt)]

    @test_inferred eltype(typeof(setcompositions_parity(3, 2))) Tuple{Tuple{SmallBitSet{UInt}, SmallBitSet{UInt}}, Bool}
    @test_inferred eltype(typeof(setcompositions_parity(SmallBitSet{UInt8}(1:5), 3, 2))) Tuple{Tuple{SmallBitSet{UInt8}, SmallBitSet{UInt8}}, Bool}
end

@testset "setcompositions" begin
    function test_compositions(ks)
        N = length(ks)
        sh = @inferred setcompositions(ks...)
        a = SmallBitSet(1:sum(ks; init = 0))
        @test @inferred(length(sh)) == factorial(big(sum(ks; init = 0))) ÷ prod(map(factorial∘big, ks); init = 1)
        @test @inferred(eltype(sh)) == NTuple{N, SmallBitSet{UInt}}
        @test all(map(length, t) == ks && (isempty(t) ? isempty(a) : (union(t...) == a)) for t in sh)
        @test allunique(sh)
    end

    function test_compositions(a::S, ks::NTuple{N,Int}) where {S <: SmallBitSet, N}
        test_compositions(ks)
        sh = @inferred setcompositions(a, ks...)
        @test @inferred(length(sh)) == factorial(big(sum(ks; init = 0))) ÷ prod(map(factorial∘big, ks); init = 1)
        @test @inferred(eltype(sh)) == NTuple{N, S}
        @test all(t isa NTuple{N, S} && map(length, t) == ks && (isempty(t) ? isempty(a) : (union(t...) == a)) for t in sh)
        @test allunique(sh)
    end

    for U in unsigned_types, (v, ks) in [(Int[], ()), (Int[], (0,)), (Int[], (0, 0)),
                (bitsize(U)-4:2:bitsize(U), (1, 1, 1)),
                (3:2:11, (5,)), (3:2:11, (2, 3)),  (3:2:11, (0, 2, 3)),  (3:2:11, (2, 0, 3)),  (3:2:11, (2, 3, 0)),
                (20:2:38, (2, 3, 2, 3)), (20:2:38, (1, 4, 0, 2, 3))]
        maximum(v; init = 0) <= bitsize(U) || continue
        a = SmallBitSet{U}(v)
        test_compositions(a, ks)
    end

    @test_throws Exception setcompositions(-1, 2)
    @test_throws Exception setcompositions(bitsize(UInt)-1, 2)
    for U in unsigned_types
        @test_throws Exception setcompositions(SmallBitSet{U}(2:2:6))
        @test_throws Exception setcompositions(SmallBitSet{U}(2:2:6), -1, 2, 2)
        @test_throws Exception setcompositions(SmallBitSet{U}(2:2:6), 3, 4)
        @test (setcompositions(SmallBitSet{U}(1:bitsize(U)), bitsize(U)-2, 2); true)
    end

    @test_inferred eltype(typeof(setcompositions(3, 2))) Tuple{SmallBitSet{UInt}, SmallBitSet{UInt}}
    @test_inferred eltype(typeof(setcompositions(SmallBitSet{UInt8}(1:5), 3, 2))) Tuple{SmallBitSet{UInt8}, SmallBitSet{UInt8}}
end

@testset "combinations(n,k)" begin
    k = 3
    @test_inferred combinations(8, k) subsets(8, k)
    s = SmallBitSet{UInt16}([2, 3, 5, 7, 8, 10, 12, 13])
    @test_inferred combinations(s, k) subsets(s, k)
    t = SmallSet{16,Int16}(s)
    @test_inferred subsets(t, k) combinations(t, k)

    sitr = subsets(s, k)
    N = 8
    T = Int8
    for C in [SmallSet{N,T}, MutableSmallSet{N,T},
            FixedVector{N,T}, MutableFixedVector{N,T}, SmallVector{N,T}, MutableSmallVector{N,T},
            PackedVector{UInt64, 5, Int8}]
        c = C(s)
        U = if C <: AbstractSmallSet
            SmallSet{N,T}
        elseif C <: Union{AbstractFixedVector, AbstractSmallVector}
            SmallVector{N,T}
        elseif C <: PackedVector
            C
        else
            error("impossible")
        end
        itr = @inferred combinations(c, k)
        @test_inferred Base.IteratorEltype(itr) Base.HasEltype()
        @test_inferred eltype(itr) U
        @test_inferred length(itr) length(sitr)
        v = collect(itr)
        @test eltype(v) == eltype(itr) && length(v) == length(sitr)
        @test Set(v) == Set(U(t) for t in sitr)
    end
end
