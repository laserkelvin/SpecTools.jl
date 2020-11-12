module SpecTools

using JLD2, Logging, ProgressMeter

export
    EnergyLevel,
    BaseLevel,
    LinearLevel,
    SymTopLevel,
    AsymTopLevel,
    Levels,
    E, g, J, K, Ka, Kc,
    Transition,
    Transitions,
    Spectrum,
    Experiment,
    StickSpectrum,
    hash_file,
    match_energy,
    match_spectrum_transitions,
    get_settings

include("base.jl")
include("constants.jl")
include("preprocess.jl")
include("utils.jl")
include("sim/Sim.jl")
include("catalog/Catalog.jl")

end
