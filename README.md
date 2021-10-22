# SpecTools

[![Build Status](https://github.com/laserkelvin/SpecTools.jl/workflows/CI/badge.svg)](https://github.com/laserkelvin/SpecTools.jl/actions)
[![Coverage](https://codecov.io/gh/laserkelvin/SpecTools.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/laserkelvin/SpecTools.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

`SpecTools` is a Julia package written to provide base abstractions for concepts in molecular spectroscopy. While rotational spectroscopy is the main usage, in principle the implementations should be general to spectroscopy at all wavelengths.

## Implementation

There are two types of representations currently implemented, which may be useful in different circumstances:

1. Concrete `Level` and `Transition` types, representing energy levels and transitions with arbitrary quantum number encoding,
2. Graph representations of levels and transitions, based on the above types.

The former is better suited for performing spectral simulations, with performant routines thanks to great packages like `Tullio.jl`, while the latter is designed specifically for doing large scale analyses and machine learning on spectroscopic graphs:

```julia
# create a bipartite graph based on a set of levels and transitions
julia> sg = BipartiteSG(levels, transitions);
# make a contiguous matrix of features for machine learning
julia> features(sg.transitions)

6×49 Matrix{Float64}:
 443.015  1148.16  2387.61  2008.85  3313.11  …  3023.19  1156.43  3436.3  2851.82
   1.0       1.0      1.0      1.0      1.0         1.0      1.0      1.0     1.0
   0.0       1.0      2.0      3.0      4.0        45.0     46.0     47.0    48.0
   0.0       0.0      0.0      0.0      0.0         0.0      0.0      0.0     0.0
   1.0       2.0      3.0      4.0      5.0        46.0     47.0     48.0    49.0
   0.0       0.0      0.0      0.0      0.0   …     0.0      0.0      0.0     0.0
```

The graph representations use `MetaGraphs.jl` to support weighted edges and vertices/nodes. For that reason, it is recommended to use `graphplot` from `GraphMakie.jl` to do visualizations of spectroscopic graphs. Similarly, we gain access to all of the graph/network analysis tools in the `JuliaGraphs` ecosystem:

```julia
julia> adjacency_matrix(bisg.graph)

99×99 SparseArrays.SparseMatrixCSC{Int64, Int64} with 192 stored entries:
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦
⠲⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠈⠳⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠈⠳⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠈⠳⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠈⠳⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠳⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠳⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
```

