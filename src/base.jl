
"""
This module defines all of the core abstraction for spectroscopy;
from energy levels to spectra.
"""

abstract type EnergyLevel{T} end

struct BaseLevel{T<:Real, R<:Integer} <: EnergyLevel{T}
    E::T
    g::R

    BaseLevel(E, g=1) = new{Float64, Int64}(E, g)
end

"""
Struct for an energy level of a linear molecule

"""
struct LinearLevel{T<:Real, R<:Integer} <: EnergyLevel{T}
    E::T
    g::R
    J::R
    F::R

    LinearLevel(E, g=1, J=1, F=0) = new{Float64, Int64}(E, g, J, F)
end

struct SymTopLevel{T<:Real, R<:Integer} <: EnergyLevel{T}
    E::T
    g::R
    J::R
    K::R
    F::R
    
    function SymTopLevel(E, g=1, J=1, K=1, F=0)
      @assert J >= K
      new{Float64, Int64}(E, g, J, K, F)
    end
end

struct AsymTopLevel{T<:Real, R<:Integer} <: EnergyLevel{T}
    E::T
    g::R
    J::R
    Ka::R
    Kc::R
    F::R

    function AsymTopLevel(E, g=1, J=1, Ka=1, Kc=1, F=0)
      @assert J >= Ka + Kc
      new{Float64, Int64}(E, g, J, Ka, Kc, F)
    end
end

# An array of energy levels
Levels = Vector{<:EnergyLevel}

# partition function terms
Q(E::Real, g::Real, T::Real) = g * exp(-E / (k_mhz * T))
Q(level::EnergyLevel, T::Real) = Q(level.E, level.g, T)
Q(levels::Levels, T::Real) = map(x -> Q(x.E, x.g, T), levels)
sumQ(levels::Levels, T::Real) = reduce(+, Q(levels, T))


struct Transition{T}
    ν::T
    ν_unc::T
    intensity::T
    lower::Union{EnergyLevel, Missing}
    upper::Union{EnergyLevel, Missing}
    Sij::T
    Aij::T

    function Transition(ν::AbstractFloat, ν_unc::AbstractFloat=0., intensity::AbstractFloat=0., lower::Union{EnergyLevel, Missing}=missing, upper::Union{EnergyLevel, Missing}=missing, Sij::AbstractFloat=0., Aij::AbstractFloat=0.)
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
