using Test, SmallCollections, SmallCombinatorics, BitIntegers

macro test_inferred(expr, good, goodtype = missing)
    msg = """

        expression:      $expr
        expected result: $good
        expected type:   $(goodtype === missing ? "type of expected result" : goodtype)
        location:        $(something(__source__.file, :none)):$(__source__.line)

        """
    quote
        let good = $good, goodtype = $goodtype,
                result = try
                    @inferred($expr)
                catch e
                    printstyled($msg; bold = true, color = :magenta)
                    rethrow(e)
                end
            if goodtype === missing
                goodtype = typeof(good)
            elseif !(goodtype isa Type)
                goodtype = typeof(goodtype)
            end
            testresult = @test isequal(result, good)
            if testresult isa Test.Pass
                testresult = @test result isa goodtype
            end
            testresult isa Test.Pass || printstyled($msg; bold = true, color = :magenta)
            result
        end
    end |> esc
end

BitIntegers.@define_integers 440

using SmallCollections: bitsize

using Random: Random, AbstractRNG, SamplerType

Random.rand(rng::AbstractRNG, ::SamplerType{Symbol}) = Symbol(rand(Char, 3)...)

unsigned_types = (UInt8, UInt64, UInt256, UInt440)

if isempty(ARGS)
    push!(ARGS, "partitions.jl", "compositions.jl", "subsets.jl", "permutations.jl")
end

foreach(include, ARGS)
