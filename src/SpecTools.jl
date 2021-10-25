module SpecTools

using Parameters
using Distributions
using Interpolations
using SHA
using PhysicalConstants
using MetaGraphs
using Tullio

# Write your package code here.

include("constants.jl")
include("types.jl")

export Level, Transition, Levels, Transitions
export e, g, ν, I, quantum_numbers, features
export match_encoding, make_linear_transitions, upper_state_qnos, lower_state_qnos

include("graphs.jl")

export BipartiteSG, ProjectedSG, offset, indices

include("lineshapes.jl")

export MultiGaussian, Gaussian, simulate_lineshape

include("signal.jl")
include("partition_function.jl")

export partition_function, read_qpart, InterpPartitionFunction

include("utils.jl")

end
