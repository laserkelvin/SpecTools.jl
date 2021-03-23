module SpecTools

using SHA, PhysicalConstants

# Write your package code here.

include("base.jl")
include("constants.jl")
include("lineshapes.jl")

export EnergyLevel,
  BaseLevel,
  LinearLevel,
  SymTopLevel,
  AsymTopLevel,
  Q,
  Experiment,
  Simulation

export Gaussian,
       Lorentzian,
       Lineshape,
       simulate_lineshape

end
