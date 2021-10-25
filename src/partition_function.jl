
# partition function terms
partition_function(E, g, T) = g * exp(-E / (k_mhz * T))

# evaluate the partition function based on levels
function partition_function(ls::Levels, temperature)
  @tullio q := partition_function(ls[i].e, ls[i].g, temperature)
end


abstract type PartitionFunction end


struct InterpPartitionFunction
  interp_obj
end

"""Constructor method for interpolated partition function objects.
Takes in a vector of temperatures, and a vector of values of the
partition function at each temperature.
"""
function InterpPartitionFunction(temperatures::Vector{T}, values::Vector{T}) where T<:Real
  interpolant = LinearInterpolation(temperatures, values)
  InterpPartitionFunction(interpolant)
end

(interp::InterpPartitionFunction)(temperature) = interp.interp_obj(temperature)

"""Read in a qpart file in the molsim format; first line is a header
that specifies the type of partition function evaluation needed, which
is parsed and used with multiple dispatch to determine how to read
the rest of the qpart file.
"""
function read_qpart(filepath)
  open(filepath) do file
    header = readline(file)
    qpart_type = occursin("interp", header) ? Symbol(IsInterp) : throw(MethodError(header, " qpart format not recognized."))
  end
  read_qpart(filepath, Val(qpart_type))
end

"""Read in an interpolated partition function
"""
function read_qpart(filepath, ::Val{:IsInterp})
  open(filepath) do file
    lines = readlines(file)
    # ignore comment lines
    data = filter(x->~startswith(lstrip(x), "#"), lines)
    # split by whitespace, then reduce into a 2xN matrix and convert to numbers
    values_matrix = mapreduce(split, hcat, data) .|> x -> parse(Float32, x)
    return InterpPartitionFunction(values_matrix[1,:], values_matrix[2,:])
  end
end
