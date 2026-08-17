@testset "ldm" begin
    for n in (1,2,6,8,14,31,48)
        s = rand(n)
        si = Int16.(rand(1:25, n))
        sv = FixedVector{n,eltype(s)}(s)
        ss = SmallSet{n,eltype(sv)}(sv)
        for k in 1:4
            @inferred ldm(s, Val(k), UInt64)
            @inferred ldm(sv, Val(k))
            @inferred ldm(ss, Val(k))
        end
        @inferred ldm(sv)
        @inferred ldm(ss)
        p = @inferred ldm(si, Val(2), NoSubsets)
        @test typeof(p) === SmallCombinatorics.LDMNode{Int16, 2, NoSubsets}
    end

    @test ldm(1:15).subsetsums === (60,60)
    @test ldm(1:15, Val(5)).subsetsums === (25, 25, 24, 23, 23)
    @test ldm(1:15, Val(5), UInt16).subsetsums === (25, 25, 24, 23, 23)
    @test ldm(1:15, Val(5), NoSubsets).subsetsums === (25, 25, 24, 23, 23)
    @test ldm(exp.(range(0,log(10),length=20)), Val(3)).value ≈ 0.1456806569078246
end
