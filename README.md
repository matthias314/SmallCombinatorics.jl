# SmallCombinatorics.jl

This package provides functions for enumerative combinatorics. It was formerly a submodule of
[SmallCollections.jl](https://github.com/matthias314/SmallCollections.jl).

The functions in `SmallCombinatorics` are often much faster than their counterparts
in other combinatorics packages. This is achieved by using the fast, non-allocating
container types offered by `SmallCollections`. These types have a maximal capacity,
which implies there are bounds on the input parameters for each function. At present,
these bounds are hard-coded. For example, the function `permutations(n)` can produce
permutations only for `n ≤ 16`.

So far, the list of implemented functions is fairly short, but hopefully it will grow.
The full documentation  for `SmallCombinatorics` is available
[here](https://matthias314.github.io/SmallCombinatorics.jl/).

## Benchmarks

The following benchmarks were done with Julia 1.11.5, Chairmarks.jl,
[Combinatorics.jl](https://github.com/JuliaMath/Combinatorics.jl) 1.0.3
and [Combinat.jl](https://github.com/jmichel7/Combinat.jl) 0.1.3.
In separate tests, GAP and Sage were 2-3 orders of magnitude slower.

### Permutations

Loop over all permutations of `1:9` and add up the image of `1` under each permutation.
The iterator returned by
[`permutations`](https://matthias314.github.io/SmallCombinatorics.jl/stable/#SmallCombinatorics.permutations)
 yields each permutation as a `SmallVector{16,Int8}`.
```julia
julia> n = 9; @b sum(p[1] for p in SmallCombinatorics.permutations($n))
1.909 ms  # 688.006 μs with @inbounds(p[1])

julia> n = 9; @b sum(p[1] for p in Combinatorics.permutations(1:$n))
14.535 s (725763 allocs: 44.297 MiB, 0.04% gc time, without a warmup)

julia> n = 9; @b sum(p[1] for p in Combinat.Permutations($n))
12.473 ms (725762 allocs: 44.297 MiB, 5.38% gc time)
```

### Combinations

Loop over all `10`-element subsets of `1:20` and add up the sum of the elements of each subset.
The iterator returned by
[`subsets`](https://matthias314.github.io/SmallCombinatorics.jl/stable/#SmallCombinatorics.subsets-Tuple{Integer,%20Integer})
yields each subset as a `SmallBitSet`.
```julia
julia> n = 20; k = 10; @b sum(sum, SmallCombinatorics.subsets($n, $k))
1.121 ms

julia> n = 20; k = 10; @b sum(sum, Combinatorics.combinations(1:$n, $k))
9.484 ms (369514 allocs: 25.373 MiB, 7.09% gc time)

julia> n = 20; k = 10; @b sum(sum, Combinat.Combinations(1:$n, $k))  # Combinat.jl
9.605 ms (369521 allocs: 25.373 MiB, 7.04% gc time)
```

### Partitions

Loop over all partitions of `40` and add up the first element of each partition.
The iterator returned by
[`partitions`](https://matthias314.github.io/SmallCombinatorics.jl/stable/#SmallCombinatorics.partitions)
 yields each permutation as a `SmallVector{64,Int8}`.
```julia
julia> n = 40; @b sum(p[1] for p in SmallCombinatorics.partitions($n))
425.938 μs  # 141.329 μs with @inbounds(p[1])

julia> n = 40; @b sum(p[1] for p in Combinatorics.integer_partitions($n))
360.947 ms (3304541 allocs: 100.824 MiB, 28.44% gc time, without a warmup)

julia> n = 40; @b sum(p[1] for p in Combinat.partitions($n))
3.761 ms (74694 allocs: 6.299 MiB)
```
