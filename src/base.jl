
"""
This module defines all of the core abstraction for spectroscopy;
from energy levels to spectra.
"""

abstract type EnergyLevel end

struct BaseLevel <: EnergyLevel
    energy::AbstractFloat
    g::Unsigned
end

struct LinearLevel <: EnergyLevel
    energy::AbstractFloat
    g::Unsigned
    J::Integer
    F::Integer
    
    # constructor method 
    function LinearLevel(energy::AbstractFloat=0., g::Unsigned=1, J::Integer=0, F::Integer=0)
      @assert g > 0
      return new(energy, g, J, F)
    end
end

struct SymTopLevel <: EnergyLevel
    energy::AbstractFloat
    g::Unsigned
    J::Integer
    K::Integer
    F::Integer
    
    # constructor method 
    function SymTopLevel(energy::AbstractFloat=0., g::Unsigned=1, J::Integer=0, K::Integer=0, F::Integer=0)
      @assert g > 0
      # angular momentum projection must always sum up to the total
      @assert J >= K
      return new(energy, g, J, K, F)
    end
end

struct AsymTopLevel <: EnergyLevel
    energy::AbstractFloat
    g::Unsigned
    J::Integer
    Ka::Integer
    Kc::Integer
    F::Integer
    
    # constructor method 
    function AsymTopLevel(energy::AbstractFloat=0., g::Unsigned=1, J::Integer=0, Ka::Integer=0, Kc::Integer=0, F::Integer=0)
      @assert g > 0
      # angular momentum projection must always sum up to the total
      @assert J >= (Ka + Kc)
      return new(energy, g, J, Ka, Kc, F)
    end
end

# An array of energy levels
Levels = Vector{EnergyLevel}

# partition function for a given state
Q(state::EnergyLevel, t::AbstractFloat) = exp(-state.E / (k * t)) * state.g

struct Transition
    ν::AbstractFloat
    ν_unc::AbstractFloat
    intensity::AbstractFloat
    lower::EnergyLevel
    upper::EnergyLevel
    Sij::AbstractFloat
    Aij::AbstractFloat

    function Transition(ν::AbstractFloat, ν_unc::AbstractFloat=0., intensity::AbstractFloat=0., lower::EnergyLevel=LinearLevel(1.), upper::EnergyLevel=LinearLevel(2.), Sij::AbstractFloat=0., Aij::AbstractFloat=0.)
      new(ν, ν_unc, intensity, lower, upper, Sij, Aij)
    end
end

Transitions = Vector{Transition}

"""
Next level of abstraction
"""

abstract type AbstractSpectrum end

struct Experiment <: AbstractSpectrum
    ν::Vector{<:AbstractFloat}
    intensity::Vector{<:AbstractFloat}
    noise::Union{Any, Vector{<:AbstractFloat}, Nothing}
end

# this is used to generate unique hashes to identify catalogs
hash_file(filepath::String) = bytes2hex(sha2_256(read(filepath)))
