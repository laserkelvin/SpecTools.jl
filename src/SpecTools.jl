module SpecTools

using Parameters, Distributions
using SHA, PhysicalConstants, MetaGraphs

# Write your package code here.

include("constants.jl")
include("types.jl")
include("graphs.jl")

export Level,
  Transition,
  Levels,
  Transitions,
  e,
  g,
  ν,
  I,
  quantum_numbers,
  features,
  match_encoding,
  make_linear_transitions,
  upper_state_qnos,
  lower_state_qnos

export BipartiteSG,
  ProjectedSG,
  offset,
  indices

using Tullio

include("lineshapes.jl")

export MultiGaussian,
  Gaussian,
  simulate_lineshape

include("signal.jl")
include("partition_function.jl")

export partition_function

include("utils.jl")

end
