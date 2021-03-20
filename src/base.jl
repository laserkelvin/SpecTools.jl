
"""
This module defines all of the core abstraction for spectroscopy;
from energy levels to spectra.
"""

abstract type EnergyLevel{T} end

struct BaseLevel{T} <: EnergyLevel{T}
    E::Float64
    g::T

    BaseLevel(E, g=1) = new{Unsigned}(E, g)
end

struct LinearLevel{T} <: EnergyLevel{T}
    E::Float64
    g::T
    J::T
    F::T

    LinearLevel(E, g=1, J=1, F=0) = new{Unsigned}(E, g, J, F)
end

struct SymTopLevel{T} <: EnergyLevel{T}
    E::Float64
    g::T
    J::T
    K::T
    F::T
    
    function SymTopLevel(E, g=1, J=1, K=1, F=0)
      @assert J >= K
      new{Unsigned}(E, g, J, K, F)
    end
end

struct AsymTopLevel{T} <: EnergyLevel{T}
    E::Float64
    g::T
    J::T
    Ka::T
    Kc::T
    F::T

    function AsymTopLevel(E, g=1, J=1, Ka=1, Kc=1, F=0)
      @assert J >= Ka + Kc
      new{Unsigned}(E, g, J, Ka, Kc, F)
    end
end

# An array of energy levels
Levels = Vector{EnergyLevel}

# partition function for a given state
Q(state::EnergyLevel, t::AbstractFloat) = exp(-state.E / (k_cm * t)) * state.g

struct Transition{T}
    ν::T
    ν_unc::T
    intensity::T
    lower::Union{EnergyLevel, Nothing}
    upper::Union{EnergyLevel, Nothing}
    Sij::T
    Aij::T

    function Transition(ν::AbstractFloat, ν_unc::AbstractFloat=0., intensity::AbstractFloat=0., lower::EnergyLevel=nothing, upper::EnergyLevel=nothing, Sij::AbstractFloat=0., Aij::AbstractFloat=0.)
      new{Float64}(ν, ν_unc, intensity, lower, upper, Sij, Aij)
    end
end

function Transition(lower::EnergyLevel, upper::EnergyLevel, intensity::AbstractFloat=0., Sij::AbstractFloat=0., Aij::AbstractFloat=0.)
  ν = upper.energy - lower.energy
  return Transition(ν, 0., intensity, lower, upper, Sij, Aij)
end

Transitions = Vector{Transition}

"""
Next level of abstraction
"""

abstract type AbstractSpectrum{T} end

struct Experiment{T} <: AbstractSpectrum{T}
    ν::Vector{<:T}
    intensity::Vector{<:T}
    noise::Union{Vector{<:T}, T, Nothing}
end

struct Simulation{T} <: AbstractSpectrum{T}
  ν::Vector{<:T}
  intensity::Vector{<:T}
  transitions::Transitions
end

# this is used to generate unique hashes to identify catalogs
hash_file(filepath::String) = bytes2hex(sha2_256(read(filepath)))
