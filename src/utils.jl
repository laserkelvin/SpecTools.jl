
"""
Convenience function to return the energy levels and
degeneracies of a linear molecule.

B is the rotational constant, D is the centrifugal
distortion term, and J_max is the maximum value
to simulate up to. Assumes that B and D are the same units.
"""
function make_linear_molecule(B, D, J_max, T=10.)
    J = 1:J_max
    E = @. B * J * (J + 1) - D * J^2 * (J + 1)^2
    G = @. 2 * (J + 1)
    levels = BaseLevel.(E, G)
    # use the boltzmann weighting as the intensity
    qrot = Q.(levels, T)
    intensity = qrot ./ sum(qrot)
    # this calculates the transition energies
    ν = diff(E)
    transitions = Transitions()
    @inbounds for i in eachindex(ν)
        trans = Transition(ν[i], 0., qrot[i+1], levels[i], levels[i+1], 0., 0.)
        push!(transitions, trans)
    end
    return levels, transitions
end


function add_noise!(I::Vector{<:Real}, σ::Real)
    target_σ = max(I) / σ
    I .+= rand(Normal(0., target_σ,), length(I))
end
