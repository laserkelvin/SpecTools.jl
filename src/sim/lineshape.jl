export
    Gaussian,
    gaussian

mutable struct Gaussian
    μ::AbstractFloat
    σ::AbstractFloat
    A::AbstractFloat
end

gaussian(x::AbstractFloat, μ::AbstractFloat, σ::AbstractFloat, A::AbstractFloat) = A * exp(-(x - μ)^2 / (2 * σ^2))

function gaussian(x::Array{<:AbstractFloat}, μ::Array{<:AbstractFloat}, σ::Array{<:AbstractFloat}, A::Array{<:AbstractFloat})
    y = zero.(x)
    for (μ_i, σ_i, A_i) in zip(μ, σ, A)
        y .+= gaussian.(x, μ_i, σ_i, A_i)
    end
    return y
end
