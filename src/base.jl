
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
    
    LinearLevel() = new()
end

LinearLevel(energy::AbstractFloat) = LinearLevel(energy, 1, 1, 1)

struct SymTopLevel <: EnergyLevel
    energy::AbstractFloat
    g::Unsigned
    J::Integer
    K::Integer
    F::Integer
    
    SymTopLevel() = new()
end

SymTopLevel(energy::AbstractFloat) = SymTopLevel(energy, 1, 1, 1, 1)

struct AsymTopLevel <: EnergyLevel
    energy::AbstractFloat
    g::Unsigned
    J::Integer
    Ka::Integer
    Kc::Integer
    F::Integer
    
    AsymTopLevel() = new()
end

AsymTopLevel(energy::AbstractFloat) = AsymTopLevel(energy, 1, 1, 1, 1, 1)

# An array of energy levels
Levels = Vector{EnergyLevel}

match_energy(level::EnergyLevel, energy::AbstractFloat) = E(level) == energy

# helper functions to grab fields
E(state::EnergyLevel) = state.energy
g(state::EnergyLevel) = state.g
J(state::EnergyLevel) = state.J
K(state::SymTopLevel) = state.K
Ka(state::AsymTopLevel) = state.Ka
Kc(state::AsymTopLevel) = state.Kc

# partition function for a given state
Q(state::EnergyLevel, t::AbstractFloat) = exp(-E(state) / (k * t)) * g(state)

struct Transition
    ν::AbstractFloat
    ν_unc::AbstractFloat
    intensity::AbstractFloat
    lower::EnergyLevel
    upper::EnergyLevel
    Sij::AbstractFloat
    Aij::AbstractFloat
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
