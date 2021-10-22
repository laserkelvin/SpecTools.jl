
# partition function terms
partition_function(E, g, T) = g * exp(-E / (k_mhz * T))

# evaluate the partition function based on levels
function partition_function(ls::Levels, temperature)
  @tullio q := partition_function(ls[i].e, ls[i].g, temperature)
end

