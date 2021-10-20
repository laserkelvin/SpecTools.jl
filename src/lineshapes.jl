
abstract type Lineshape{T} end

"""Compute the amplitude of a Gaussian lineshape, given a grid value
x, and parameters A, μ, σ as the amplitude, center, and width
respectively.
"""
gaussian(x, (A, μ, σ)) = A / (σ*(2π)^0.5) * exp(-(x - μ)^2 / (2σ^2))
gaussian(x, A, μ, σ) = A / (σ*(2π)^0.5) * exp(-(x - μ)^2 / (2σ^2))

"""A struct of arrays for holding multiple Gaussian lineshapes.
"""
struct MultiGaussian{T} <: Lineshape{T}
  A::Vector{T}
  μ::Vector{T}
  σ::Vector{T}
end

struct Gaussian{T} <: Lineshape{T}
  A::T
  μ::T
  σ::T
end

"""Constructor method for multiple Gaussians, using a matrix
of parameters.
"""
MultiGaussian(θ) = MultiGaussian(θ[1,1:end], θ[2,1:end], θ[3,1:end])

# functor definitions
function (g::MultiGaussian)(x::Real)
  y = zero(x)
  @inbounds for i ∈ eachindex(g.A, g.μ, g.σ)
    y += gaussian(x, (g.A[i], g.μ[i], g.σ[i]))
  end
end
(g::Gaussian)(x) = gaussian(x, (g.A, g.μ, g.σ))

# now for more specialized, performant variants

"""For a vector of x (length M), we can use einsum notation
to minimize allocations and get a significant boost in speed.

The notation is such that y[m] (same length as x) requires
summing over N, the number of Gaussians.

This algorithm is still O(MN), and so improvements could be
made with fast Gauss transforms.
"""
function (g::MultiGaussian)(x::Vector{T}) where T<:Real
  A, μ, σ = g.A, g.μ, g.σ
  @tullio y[m] := A[n] / (2*π*σ[n]^2)^0.5 * exp(-(x[m] - μ[n])^2 / (2*σ[n]^2))
end

"""Given a grid value and a generic lineshape struct, this function
provides a general interface to simulate any lineshape that is
a concrete type of Lineshape.
"""
simulate_lineshape(x, l::Lineshape) = l(x)

