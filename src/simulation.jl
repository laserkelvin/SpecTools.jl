
# calculate the energy level population
population(g, E, T) = g * exp(-E / (k_B * T))

function τ(aij, n, g_up, E_up, T, ν, δv, Q)
    value = (aij * c_0^3 * (n * 100^2)) * population(g_up, E_up, T)
    value *= (exp((h * ν * 1e6) / (k_B * T)) - 1)
    value /= (8 * π * (ν * 1e6)^3 * δv * 1000 * Q)
end

# this is the expression for brightness temperature
function Tb(ν, T)
    hν = h * ν * 1e6
    J_T = hν * (exp((hν) / (k_B * T)) - 1)^(-1)
    return J_T
end

# This is the excess brightness temperature, according to Turner (1991)
delta_Tb(ν, τ, T_ul, T_bg) = (Tb(ν, T_ul) - Tb(ν, T_bg)) * (1 - exp(-τ))

beam_size(ν, dish_size) = (206265. * 1.22 * (c_0 / (ν * 1e6))) / dish_size

# perform an in-place correction to the intensity, 
# written to minimize allocations
function beam_correction!(ν, I, source_size, dish_size)
    # since we disable bounds checking, we should be sure
    # we don't access out of bounds
    @assert sizeof(ν) == sizeof(I) "Frequency and intensity array sizes don't match!"
    # disable bounds checking for performance
    @inbounds for idx in eachindex(I)
        I[idx] *= (source_size^2 / (beam_size(ν[idx], dish_size)^2 + source_size^2))
    end
end