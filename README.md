# SmallCombinatorics.jl

This package provides functions for enumerative combinatorics. It was formerly a submodule of
[SmallCollections.jl](https://github.com/matthias314/SmallCollections.jl).

The functions in `SmallCombinatorics` are often much faster than their counterparts
in other combinatorics packages. This is achieved by using the fast, non-allocating
container types offered by `SmallCollections`. These types have a maximal capacity,
which implies there are bounds on the input parameters of the functions. At present,
most bounds are hard-coded. For example, the function `permutations(n::Integer)` can produce
permutations only for `n ≤ 16`. In contrast, the vector version `permutations(v::SmallVector)`
has no restriction on the length of the
[`SmallVector`](https://matthias314.github.io/SmallCollections.jl/v0.5.0/smallvector/#SmallCollections.SmallVector).

So far, the list of implemented functions is fairly short, but hopefully it will grow.
The full documentation  for `SmallCombinatorics` is available
[here](https://matthias314.github.io/SmallCombinatorics.jl/).

## Benchmarks

The following benchmarks were done with Julia 1.11.6, Chairmarks.jl, SmallCombinatorics.jl v0.1.1,
[Combinatorics.jl](https://github.com/JuliaMath/Combinatorics.jl) 1.0.3
and [Combinat.jl](https://github.com/jmichel7/Combinat.jl) 0.1.4.
In separate tests, GAP and Sage were 2-3 orders of magnitude slower.

### Permutations

Loop over all permutations of `1:9` and add up the image of `1` under each permutation.
The iterator returned by
[`SmallCombinatorics.permutations`](https://matthias314.github.io/SmallCombinatorics.jl/stable/#SmallCombinatorics.permutations)
 yields each permutation as a `SmallVector{16,Int8}`.
```julia
julia> n = 9; @b sum(@inbounds(p[1]) for p in SmallCombinatorics.permutations($n))
683.524 μs

julia> n = 9; @b sum(@inbounds(p[1]) for p in Combinatorics.permutations(1:$n))
14.171 s (725763 allocs: 44.297 MiB, 0.40% gc time, without a warmup)

julia> n = 9; @b sum(@inbounds(p[1]) for p in Combinat.Permutations($n))
13.309 ms (725762 allocs: 44.297 MiB, 11.81% gc time)
```

### Combinations

Loop over all `10`-element subsets of `1:20` and add up the sum of the elements of each subset.
The iterator returned by
[`SmallCombinatorics.combinations`](https://matthias314.github.io/SmallCombinatorics.jl/stable/#SmallCombinatorics.combinations-Tuple{Integer,%20Integer})
yields each subset as a [`SmallBitSet`](https://matthias314.github.io/SmallCollections.jl/stable/smallbitset/#SmallCollections.SmallBitSet).
```julia
julia> n = 20; k = 10; @b sum(first(c) for c in SmallCombinatorics.combinations($n, $k))
378.624 μs

julia> n = 20; k = 10; @b sum(first(c) for c in Combinatorics.combinations(1:$n, $k))
9.064 ms (369514 allocs: 25.373 MiB, 6.34% gc time)

julia> n = 20; k = 10; @b sum(first(c) for c in Combinat.Combinations(1:$n, $k))
7.768 ms (184765 allocs: 19.735 MiB, 2.41% gc time)
```

### Partitions

Loop over all partitions of `40` and add up the first element of each partition.
The iterator returned by
[`SmallCombinatorics.partitions`](https://matthias314.github.io/SmallCombinatorics.jl/stable/#SmallCombinatorics.partitions)
 yields each permutation as a `SmallVector{64,Int8}`.
```julia
julia> n = 40; @b sum(@inbounds(p[1]) for p in SmallCombinatorics.partitions($n))
141.197 μs

julia> n = 40; @b sum(@inbounds(p[1]) for p in Combinatorics.integer_partitions($n))
316.759 ms (3304541 allocs: 100.823 MiB, 24.99% gc time, without a warmup)

julia> n = 40; @b sum(@inbounds(p[1]) for p in Combinat.Partitions($n))
2.596 ms (37340 allocs: 4.366 MiB)
```
