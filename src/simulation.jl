
# Expression for the optical depth of a transition
function τ(aij, n, g_up, E_up, T, ν, δv, qsum)
    value = (aij * c_0^3 * (n * 100^2)) * Q(E_up, g_up, T)
    value *= (exp((h * ν * 1e6) / (k_B * T)) - 1)
    value /= (8 * π * (ν * 1e6)^3 * δv * 1000 * qsum)
end

function τ(t::Transition, n, T, δv, qsum)
    aij, g_up, E_up, ν = t.Aij, t.upper.g, t.upper.E, t.ν
    return τ(aij, n, g_up, E_up, T, ν, δv, qsum)
end

# this is the expression for brightness temperature
function Tb(ν, T)
    hν = h * ν * 1e6
    J_T = hν * (exp((hν) / (k_B * T)) - 1)^(-1)
    return J_T
end

# This is the excess brightness temperature, according to Turner (1991)
ΔTb(ν, τ, T_ul, T_bg) = (Tb(ν, T_ul) - Tb(ν, T_bg)) * (1 - exp(-τ))

# compute the beam size, given a frequency and the size of the dish
beam_size(ν, dish_size) = (206265. * 1.22 * (c_0 / (ν * 1e6))) / dish_size

# perform an in-place correction to the intensity, 
# written to minimize allocations
function beam_correction!(ν, I, source_size, dish_size)
    # disable bounds checking for performance
    @inbounds for i in eachindex(I, ν)
        I[i] *= (source_size^2 / (beam_size(ν[i], dish_size)^2 + source_size^2))
    end
end
