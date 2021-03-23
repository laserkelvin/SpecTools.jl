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
  Q,
  Experiment,
  Simulation

include("lineshapes.jl")

export Gaussian,
       Lorentzian,
       Lineshape,
       simulate_lineshape

include("signal.jl")
include("simulation.jl")

end
