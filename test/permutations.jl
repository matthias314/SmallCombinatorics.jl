using SmallCombinatorics: PermN, PermEltype

@testset "permutations(n)" begin
    @test_throws Exception permutations(-1)
    @test_throws Exception permutations(PermN+1)

    for n in [0, 1, 2, 7]
        p = @inferred permutations(n)
        @test_inferred eltype(p) SmallVector{PermN,PermEltype}
        @test_inferred length(p) factorial(n)
        c = collect(p)
        @test length(c) == length(p) && eltype(c) == eltype(p)
        @test allunique(c) && all(w -> Set(w) == Set(1:n), c)
    end
end

VERSION >= v"1.11" && @testset "permutations(v)" begin
    for V in [FixedVector, MutableFixedVector], N in [1, 2, 7], T in [Int8, UInt32, Char, Symbol]
        local v
        while true
            v = rand(V{N,T})
            allunique(v) && break
        end
        p = @inferred permutations(v)
        @test_inferred eltype(p) SmallVector{N,T}
        @test_inferred length(p) factorial(N)
        c = collect(p)
        @test length(c) == length(p) && eltype(c) == eltype(p)
        @test allunique(c) && all(w -> Set(w) == Set(v), c)

        q = @inferred permutations_parity_transposition(v)
        @test_inferred Base.IteratorEltype(q) Base.HasEltype()
        @test_inferred eltype(q) Tuple{SmallVector{N,T}, Bool, Tuple{Int,Int}}
        @test_inferred length(q) factorial(N)
    end

    for V in [SmallVector, MutableSmallVector], N in [1, 2, 7], T in [Int8, UInt32, Char, Symbol]
        local v
        while true
            v = rand(V{N,T})
            allunique(v) && break
        end
        p = @inferred permutations(v)
        @test_inferred eltype(p) SmallVector{N,T}
        @test_inferred length(p) factorial(length(v))
        c = collect(p)
        @test length(c) == length(p) && eltype(c) == eltype(p)
        @test allunique(c) && all(w -> Set(w) == Set(v), c)

        q = @inferred permutations_parity_transposition(v)
        @test_inferred Base.IteratorEltype(q) Base.HasEltype()
        @test_inferred eltype(q) Tuple{SmallVector{N,T}, Bool, Tuple{Int,Int}}
        @test_inferred length(q) factorial(length(v))
    end

    local v
    T = Int16
    while true
        v = rand(PackedVector{UInt32,5,T})
        allunique(v) && break
    end
    N = capacity(v)
    p = @inferred permutations(v)
    @test_inferred eltype(p) SmallVector{N,T}
    @test_inferred length(p) factorial(length(v))
    c = collect(p)
    @test length(c) == length(p) && eltype(c) == eltype(p)
    @test allunique(c) && all(w -> Set(w) == Set(v), c)

    q = @inferred permutations_parity_transposition(v)
    @test_inferred Base.IteratorEltype(q) Base.HasEltype()
    @test_inferred eltype(q) Tuple{SmallVector{N,T}, Bool, Tuple{Int,Int}}
    @test_inferred length(q) factorial(length(v))

    v = SmallVector{typemax(PermEltype)+1,Int16}(0:typemax(PermEltype))
    @test_throws Exception permutations(v)
end

VERSION >= v"1.11" && @testset "multi_permutations" begin
    for V in [SmallVector, MutableSmallVector, FixedVector, MutableFixedVector], N in [1, 2, 3], T in [Bool, Int8, Char]
        T == Bool && N > 2 && continue
        local w
        while true
            w = rand(V{N,T})
            allunique(w) && break
        end
        v = V{2*N,T}([w..., w...])
        p = @inferred multiset_permutations(v)
        d = unique(permutations(v))
        @test_inferred eltype(p) SmallVector{2*N,T}
        @test_inferred length(p) length(d)
        c = collect(p)
        @test length(c) == length(d) && eltype(c) == eltype(d)
        @test allunique(c) && Set(c) == Set(d)
    end
end
