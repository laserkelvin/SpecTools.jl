module SpecTools

using SHA, PhysicalConstants

# Write your package code here.

include("base.jl")
include("constants.jl")

export EnergyLevel,
  BaseLevel,
  LinearLevel,
  SymTopLevel,
  AsymTopLevel,
  Levels,
  Q,
  sumQ,
  Experiment,
  Simulation,
  Transition,
  Transitions

include("lineshapes.jl")

export MultiGaussian,
  Gaussian,
  simulate_lineshape

include("signal.jl")
include("simulation.jl")

include("utils.jl")

export make_linear_molecule

end
