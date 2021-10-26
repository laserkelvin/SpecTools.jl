
const labels = ["A1", "A2", "E", "B1", "B2", "A"]
const label_map = Dict(key => "$num" for (num, key) ∈ enumerate(labels))
# add the extra character
label_map["-"] = ""
const column_types = [Int64, Int64, Int64, Float64, Int64, Float64, Float64]

"""Function to parse in a CSV formatted line list output from PGopher.
Kwargs specify what to skip and where to read the header in: by default,
the options expect to skip the first line, read the header from the second,
startin parsing from the third, and skips the last line.
"""
function read_pgopher_linelist(filepath; header=2, skipto=3, footerskip=1)
    csv_data = CSV.File(filepath, header=header, skipto=skipto, footerskip=footerskip)
    ν = csv_data[:Position]
    A = csv_data[:A]
    encodings = clean_linelist_encoding(csv_data)
    transitions = StructArray([Transition(ν=ν[i], I=A[i], encoding=encodings[:,i]) for i ∈ eachindex(ν, A)])
    return transitions
end

"""Parse through an energy level list output of PGopher.
"""
function read_pgopher_levels(filepath)
    open(filepath) do file
        lines = readlines(file)[3:end-3]
        return StructArray(parse_level_row.(lines))
    end
end

"""Not intended for public consumption, but encapsulating the work done
on a row actually helps the eventual `StructArray` creation infer
the correct types, which should give downstream performance improvements
(instead of just `Vector{Any}`).
"""
function parse_level_row(row)
    split_line = replace_sym_labels(row) |> split
    # for remaining columns, parse as integers
    remaining = abs(length(column_types) - length(split_line))
    expected_types = vcat(copy(column_types), [Int64 for i in 1:remaining])
    data = [parse(expec_type, value) for (expec_type, value) in zip(expected_types, split_line)]
    Level(e=data[4], g=data[5], encoding=[val for val in data[8:end]])
end

function replace_sym_labels(string_encoding)
    for (key, value) ∈ label_map
        string_encoding = replace(string_encoding, key=>value)
    end
    string_encoding
end

"""Function to clean up the linelist quantum number encoding. The
first action replaces symmetry labels with a numeric mapping, followed
by splitting each row into a vector of strings. The result is combined
to yield a `M x N` matrix for `M` quantum numbers, and `N` entries.

Finally, the matrix values are converted into integers.
"""
function clean_linelist_encoding(csv_data)
    mapreduce(x->replace_sym_labels(x) |> split, hcat, csv_data[:Label]) |> x -> parse.(Int, x)
end
