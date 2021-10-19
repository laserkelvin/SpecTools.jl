
abstract type Lineshape{T} end

"""Compute the amplitude of a Gaussian lineshape, given a grid value
x, and parameters A, μ, σ as the amplitude, center, and width
respectively.
"""
gaussian(x, (A, μ, σ)) = A / (σ*(2π)^0.5) * exp(-(x - μ)^2 / (2σ^2))

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
(g::MultiGaussian)(x) = mapreduce(gaussian, +, Ref(x), (g.A, g.μ, g.σ))
(g::Gaussian)(x) = gaussian(x, (g.A, g.μ, g.σ))

"""Given a grid value and a generic lineshape struct, this function
provides a general interface to simulate any lineshape that is
a concrete type of Lineshape.
"""
simulate_lineshape(x, l::Lineshape) = l(x)

function decompose_multigaussians(x, g::MultiGaussian)
  map(gaussian, Ref(x), (g.A, g.μ, g.σ))
end

