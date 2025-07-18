```@meta
DocTestSetup = :(using SmallCollections, SmallCombinatorics)
```

# SmallCombinatorics.jl

```@docs
SmallCombinatorics
```
# Partitions
```@docs
partitions
```

# Compositions
```@docs
compositions
compositions_cumsum
weakcompositions
weakcompositions_cumsum
```

## Subsets and set compositions

When used with a `SmallBitSet` as first argument, the following functions internally use
the function `SmallCollections.pdep`.
As discussed in the docstring for `pdep`, performance is much better if the processor supports the BMI2 instruction set.
The same applies to `setcompositions` with more than two parts, even if the first argument is not a `SmallBitSet`.

```@docs
subsets(::Integer)
subsets(::Integer, ::Integer)
setcompositions
setcompositions_parity(::Vararg{Integer})
setcomposition_parity
```

## Combinations

```@docs
combinations(::Integer, ::Integer)
```

## Permutations

```@docs
permutations
permutations_parity_transposition
```
