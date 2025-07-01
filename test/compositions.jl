@testset "compositions" begin
    using SmallCombinatorics: CompN, CompEltype

    for n in [-1, 0, 8], k in [-1, 0, 4, 8, 9, CompN+1]
        if n < 0 || k < 0 || k > CompN
            @test_throws Exception weakcompositions(n, k)
            @test_throws Exception compositions(n, k)
            continue
        end
        @test_inferred eltype(typeof(weakcompositions(n, k))) SmallVector{CompN,CompEltype}
        v = @inferred collect(weakcompositions(n, k))
        @test_inferred eltype(typeof(compositions(n, k))) SmallVector{CompN,CompEltype}
        w = @inferred collect(compositions(n, k))
        if k == n == 0
            @test isempty(only(v))
        else
            @test length(v) == length(weakcompositions(n, k)) == binomial(n+k-1, k-1)
            @test all(c -> all(>=(0), c) && length(c) == k, v)
            @test all(c -> sum(c) == n, v)
            @test allunique(v)

            @test length(w) == length(compositions(n, k)) == (k <= n ? binomial(n-1, k-1) : 0)
            @test all(c -> all(>(0), c) && length(c) == k, w)
            @test all(c -> sum(c) == n, w)
            @test allunique(w)
        end
    end
end
