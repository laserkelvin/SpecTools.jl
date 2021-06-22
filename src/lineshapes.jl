
abstract type Lineshape{T} end


struct Gaussian{T} <: Lineshape{T}
  μ::T
  σ::T
  A::T
end

# constructor method from an array
function Gaussian(p::Array{<:Float64,2})
  return Gaussian.(p[:,1], p[:,2], p[:,3])
end

# general expression for a Gaussian
function simulate_lineshape(x::AbstractFloat, g::Gaussian)
    abs(x - g.μ) > 20 * g.σ && return 0.
    return (g.A / (2 * π * g.σ^2)) * exp(-(x - g.μ)^2 / (2 * g.σ^2))
end

struct Lorentzian{T} <: Lineshape{T}
  μ::T
  fwhm::T
  A::T
end

# constructor method from an array
function Lorentzian(p::Array{<:Float64,2})
  return Lorentzian.(p[:,1], p[:,2], p[:,3])
end

function simulate_lineshape(x::AbstractFloat, l::Lorentzian)
    return (l.A / π) * (l.fwhm / (x - l.μ)^2 + l.fwhm^2)
end

"""
  Generalized functions for the Lineshape supertype
"""

function (lineshape::Lineshape)(x)
  simulate_lineshape(x, lineshape)
end


function (lineshapes::Vector{<:Lineshape})(x)
  y = 0.
  for lineshape in lineshapes
    y += lineshape(x)
  end
  return y
end


Gaussian(t::Transition, width::Real) = Gaussian(t.ν, width, t.intensity)

Lorentzian(t::Transition, width::Real) = Lorentzian(t.ν, width, t.intensity)
