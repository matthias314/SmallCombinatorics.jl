@testset "ldm" begin
    for n in (1,2,6,8,48,120,1243)
        s = rand(n)
        for k in (2,3,4)
            ldm(s, Val(k))
        end
    end

    @test ldm(1:11).subsetsums === (33, 33)
    @test ldm(1:11, Val(4)).subsetsums === (18, 17, 16, 15)
    @test ldm(1:11, Val(4), UInt16).subsetsums === (18, 17, 16, 15)
    @test ldm(1:11, Val(4), NoSubsets).subsetsums === (18, 17, 16, 15)
    @test ldm(exp.(range(0,log(10),length=20)), Val(3)).value ≈ 0.1456806569078246
end
