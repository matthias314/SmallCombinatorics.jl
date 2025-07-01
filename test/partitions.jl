@testset "partitions" begin
    using SmallCombinatorics: PartN, PartEltype

    for n in [-1, 0, 8, PartN+1]
        if n < 0 || n > PartN
            @test_throws Exception partitions(n)
            continue
        end
        @test_inferred eltype(typeof(partitions(n))) SmallVector{PartN,PartEltype}
        v = @inferred collect(partitions(n))
        @test all(p -> issorted(p; rev = true) && all(>(0), p), v)
        @test all(p -> sum(p; init = 0) == n, v)
        @test allunique(v)
        n == 0 && @test length(v) == 1
        n == 8 && @test length(v) == 22
    end

    for n in [-1, 0, 20, typemax(PartEltype)+1], k in [-1, 0, 1, 10, 20, PartN+1]
        if !(0 <= n <= typemax(PartEltype) && 0 <= k <= PartN)
            @test_throws Exception partitions(n, k)
            continue
        end
        @test_inferred eltype(typeof(partitions(n, k))) SmallVector{PartN,PartEltype}
        v = @inferred collect(partitions(n, k))
        @test all(p -> issorted(p; rev = true) && all(>(0), p), v)
        @test all(p -> length(p) == k && sum(p; init = 0) == n, v)
        @test allunique(v)
        p20 = Dict(0 => 0, 1 => 1, 10 => 42, 20 => 1)
        n == 20 && @test length(v) == p20[k]
    end
end
