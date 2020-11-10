
export
    Q,
    calculate_τ,
    vlsr,
    calculate_beam_dilution,
    apply_beam_corrections!

include("catalog/base.jl")

"""
Expressions for the partition function. Basically, degeneracy multiplied by
Boltzmann factor at a given temperature.
"""
Q(g::Integer, E::AbstractFloat, T::AbstractFloat) = g * exp(-(E / (c_0.val / 1e4)) / (k_cm * T))

# vectorized calculation of partition function for a vector of levels
Q(levels::Vector{<:EnergyLevel}, T::AbstractFloat) = sum(Q.(g.(levels), E.(levels), T))

function calculate_τ(transition::Transition, molecule::Molecule, dV::AbstractFloat)
    # grab some fields
    g_up, e_up = g(transition.upper), E(transition.upper)
    ν, Aij = transition.ν, transition.Aij
    ncol, T, Q = molecule.column_density, molecule.temperature, molecule.partition_function
    # actually do the calculation now
    τ = Aij * c_0.val^3. * (ncol * 100^2) * g_up * exp(-e_up / T)
    τ *= exp(h.val * ν * 1e6 / (k_B * T)) - 1
    τ /= (8 * π * (ν * 1e6)^3 * dV * 1000 * Q)
    return τ
end

"""
Calculate the Doppler shifted frequency, given a rest frequency
ν₀ and a radial velocity in km/s
"""
function vlsr(ν₀::AbstractFloat, velocity::AbstractFloat)
    return ν₀ - (velocity * ν₀ / (c_0.val / 1e3))
end

"""
Calculate the beam dilution factor at a given frequency, source and dish sizes.
"""
function calculate_beam_dilution(frequency::AbstractFloat, source_size::AbstractFloat, dish_size::AbstractFloat)
    beam_size = 206265 * 1.22 * (c_0.val / frequency * 1e6) / dish_size
    dilution_factor = source_size^2 / (beam_size^2 + source_size^2)
    return dilution_factor
end

"""
Vectorized inplace correction of the simulated flux, given source and dish sizes.
"""
function apply_beam_corrections!(s::Spectrum, source_size::AbstractFloat, dish_size::AbstractFloat)
    s.intensity .*= calculate_beam_dilution.(s.frequency, source_size, dish_size)
    return
end