
module Sim

using ..SpecTools

export
    Gaussian,
    gaussian

include("lineshape.jl")

export
    Q,
    calculate_τ,
    vlsr,
    calculate_beam_dilution,
    apply_beam_corrections!,
    doppler_shift_ν

include("compute.jl")

end