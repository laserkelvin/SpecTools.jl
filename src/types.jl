
"""This file contains basic abstractions for molecular spectroscopy
in Julia, based around the two types `Level` and `Transition`
corresponding to energy levels, and transitions between energy
levels respectively.
"""


"""Concrete type for an energy level. The three fields
correspond to the energy of the level, its degeneracy,
and a vector of integers encoding its quantum numbering.

For performance reasons, the encoding is stored as a
vector as opposed to a dictionary mapping with human
readable names. The `quantum_numbers` function provides
an interface to obtain this.
"""
@with_kw struct Level{T,U,S}
	e::T = 1f0
	g::U = 1
	encoding::Vector{S} = [0]
end

"""Concrete type for a transition object. The three fields
correspond to the frequency of the transition, the intensity
of the transition (such as the intrinsic linestrength, the
Einstein A coefficient), and the quantum number encoding
for the transition.

While not explicitly enforced, the encoding is expected to
contain an even number of elements corresponding to lower
and upper state quantum numbers.
"""
@with_kw struct Transition{T,U,S}
	ν::T = 1f0
	I::U = 1f0
	encoding::Vector{S} = [0,1]
    lower::Union{Nothing, Level} = nothing
    upper::Union{Nothing, Level} = nothing
end

# Aliases for collections of either types
Levels = Vector{Level}
Transitions = Vector{Transition}

e(l::Level) = l.e
g(l::Level) = l.g

ν(t::Transition) = t.ν
I(t::Transition) = t.I

# Generate a mapping between string quantum numbers and their value
function quantum_numbers(l::Level, names...)
	return Dict(Symbol(name) => value for (name, value) in zip(names, l.encoding))
end

function quantum_numbers(t::Transition, names...)
	encoding = Dict()
	for (name, l, u) in zip(names, lower_state_qnos(t), upper_state_qnos(t))
		encoding = merge(encoding, Dict(Symbol("$name", "_l") => l, Symbol("$name", "_u") => u))
	end
	return encoding
end

"""Function that converts the fields in a `Level` object into
a tuple of a single, promoted type.
"""
features(l::Level) = promote([l.e, l.g, l.encoding...]...)
features(t::Transition) = promote([t.ν, t.I, t.encoding...]...)

"""Function that converts a collection of spectroscopic types
into a contiguous matrix of promoted types, with the intention
of doing some kind of computation on them.
"""
features(ts) = mapreduce(x -> collect(features(x)), hcat, ts)

# function to compute the line intensity from Aij and partition function
intensity(t::Transition, z::Real) = t.I / z

"""Function to randomly generate linear molecule energy levels.
This is generally for testing purposes than for actual use.
"""
function make_linear_levels(n)
	e = rand(Uniform(0f0, 1f5), n)
	sort!(e)
	j = [[i, 0] for i ∈ 0:n-1]
	g = [2 * i + 1 for i ∈ 0:n-1]
	levels = Levels()
	sizehint!(levels, n)
	for i ∈ eachindex(e,j,g)
		push!(levels, Level(g=g[i], e=e[i], encoding=j[i]))
	end
	levels
end

"""Find the index of an spectroscopic object from a vector
of objects that matches the encoding provided. In other words,
match the quantum numbers.
"""
function match_encoding(obj, enc)
	@inbounds for i ∈ eachindex(obj)
		if obj[i].encoding == enc
			return i
		end
	end
end

"""Function that processes an iterable of energy levels, and
creates a vector of transition objects from them, assuming that
transitions between each contiguous energy level.

Mostly for testing purposes.

A quick way to do this would be to pipe the output of
`make_linear_levels` into this function.
"""
function make_linear_transitions(levels)
    transitions = Transitions()
	sizehint!(transitions, length(levels))
	for i ∈ 1:length(levels)-1
		l, u = levels[i], levels[i+1]
		ν = abs(e(l) - e(u))
		encoding = vcat(l.encoding, u.encoding)
		push!(transitions, Transition(ν=ν, encoding=encoding, lower=l, upper=u))
	end
	return transitions
end

function upper_state_qnos(t::Transition)
	offset = length(t.encoding) ÷ 2
	return t.encoding[offset+1:length(t.encoding)]
end

function lower_state_qnos(t::Transition)
	offset = length(t.encoding) ÷ 2
	return t.encoding[1:offset]
end

