
export
    Q

"""
Computation
"""


"""
Expressions for the partition function. Basically, degeneracy multiplied by
Boltzmann factor at a given temperature.
"""
Q(g::Integer, E::AbstractFloat, T::AbstractFloat) = g * exp(-(E / (c_0.val / 1e4)) / (k_cm * T))

# vectorized calculation of partition function for a vector of levels
Q(levels::Vector{<:EnergyLevel}, T::AbstractFloat) = sum(Q.(g.(levels), E.(levels), T))