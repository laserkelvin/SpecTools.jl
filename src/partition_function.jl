
# partition function terms
Q(E, g, T) = g * exp(-E / (k_mhz * T))
#Q(level::EnergyLevel, T::Real) = Q(level.E, level.g, T)
#Q(levels::Levels, T::Real) = map(x -> Q(x.E, x.g, T), levels)
#sumQ(levels::Levels, T::Real) = reduce(+, Q(levels, T))
#
