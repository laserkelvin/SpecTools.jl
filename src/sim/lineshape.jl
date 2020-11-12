
abstract type Lineshape end

mutable struct Gaussian <: Lineshape
    μ::AbstractFloat
    σ::AbstractFloat
    A::AbstractFloat
end

gaussian(x::AbstractFloat, μ::AbstractFloat, σ::AbstractFloat, A::AbstractFloat) = A * exp(-(x - μ)^2 / (2 * σ^2))

# vectorized version for multiple peaks
function gaussian(x::Vector{<:AbstractFloat}, μ::Vector{<:AbstractFloat}, σ::AbstractFloat, A::Vector{<:AbstractFloat})
    y = zero.(x)
    for (μ_i, A_i) in zip(μ, A)
        y .+= gaussian.(x, μ_i, σ, A_i)
    end
    return y
end
