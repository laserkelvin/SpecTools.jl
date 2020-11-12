
using FileIO, ProgressMeter

function read_pgopher_levels(filepath::String)
    # read in the whitespace delimited file, skip first two lines
    df = DataFrame(readdlm(filepath, skipstart=2))
    # get rid of bad lines
    mask = isa.(df[:,1], Number)
    df = df[mask,:]
    # types never inferred correctly
    df[!,[4, 6]] = convert.(Float64, df[:,[4,6]])
    df[!,[1, 3, 5, 8, 9, 10, 12]] = convert.(Integer, df[:,[1, 3, 5, 8, 9, 10, 12]])
    # extracts the bare minimum information, and rename
    df = select(df, :x4 => :energy, :x5 => :g)
    return df
end

"""
Function that calls the command line PGopher to calculate the
energy levels.
"""
function run_pgopher_levels(filepath::String, nproc::Integer=4)
    contents = read(`pgo --np $nproc --energies $filepath`, String)
    lines = split(contents, "\n")[3:end-3]
    energies, g = Vector{Float64}(), Vector{Int64}()
    for (index, line) in enumerate(lines)
        contents = split(line)
        push!(energies, parse(Float64, contents[4]))
        push!(g, parse(Int64, contents[5]))
    end
    return DataFrame(energy=energies, g=g)
end

function run_pgopher_catalog(filepath::String, nproc::Integer=4)
    contents = read(`pgo --np $nproc $filepath`, String)
    lines = split(contents, "\n")[2:end-2]
    results = []
    for line in lines
        if occursin(r"[A-Z]\d\sX", line)
            split_line = split(line)
            indices = [10, 11, 12, 13, 14, 15]
            qnos = [20, 21, 22, 24, 26, 27, 28, 30]
            # convert strings to floats
            values = parse.(Float64, split_line[indices])
            push!(results, values)
        end
    end
    dataframe = DataFrame(transpose(reduce(hcat, results)), ["ν", "intensity", "e_up", "e_low", "Sij", "Aij"])
    return dataframe
end

function read_pgopher_catalog(filepath::String)
    open(filepath) do file
        results = []
        for line in eachline(file)
            # use regex to find all lines that have transitions data
            if occursin(r"[A-Z]\d\sX", line)
                split_line = split(line)
                indices = [10, 11, 12, 13, 14, 15]
                # qnos = [20, 21, 22, 24, 26, 27, 28, 30]
                # convert strings to floats
                values = parse.(Float64, split_line[indices])
                push!(results, values)
            end
        end
        dataframe = DataFrame(transpose(reduce(hcat, results)), ["ν", "intensity", "e_up", "e_low", "Sij", "Aij"])
        return dataframe
    end
end

"""
Run command line PGopher to extract energy levels and the transition
catalog. This assumes that the correct settings are provided in the .pgo
file, for example that EinsteinA is desired, the temperature to simulate
at, and all the other knick knacks.
"""
function run_pgopher(filepath::String, nproc::Integer=4)
    levels = run_pgopher_levels(filepath, nproc)
    # levels = select(levels, :energy => :e_low, :g)
    catalog = run_pgopher_catalog(filepath, nproc)
    # merge the data into a single dataframe
    return levels, catalog
end

function merge_pgopher_results(levels_df::DataFrame, catalog_df::DataFrame)
    upper = select(levels_df, :energy => :e_up, :g => :g_up)
    lower = select(levels_df, :energy => :e_low, :g => :g_low)
    # merge dataframes to give every transition the necessary information
    combined = innerjoin(catalog_df, lower, on = :e_low)
    combined = innerjoin(combined, upper, on = :e_up)
    return combined
end

function dataframes_to_objects(levels_df::DataFrame, combined_df::DataFrame)
    transitions = Transitions()
    # sort the levels dataframe in ascending energy
    sort!(levels_df, [:energy])
    # vectorized creation of EnergyLevels
    levels = BaseLevel.(levels_df.energy, levels_df.g)
    # now that all the level objects are created, we can link transitions
    # to lookups for each level
    p = Progress(size(combined_df)[1], 1, "Generating Transition objects...", 50)
    for row in eachrow(combined_df)
        next!(p)
        # match the levels up
        low_idx, up_idx = searchsortedfirst(levels_df.energy, row.e_low), searchsortedfirst(levels_df.energy, row.e_up)
        push!(
            transitions,
            Transition(
                row.ν,
                0.,
                row.intensity,
                levels[low_idx], levels[up_idx],
                row.Sij,
                row.Aij
            )
            )
    end
    return levels, transitions
end

struct PGopher <: BaseCatalog
    transitions::Vector{<:Transition}
    levels::Vector{<:EnergyLevel}
    name::String
    hash::String
end

"""
Construct a `PGopher` object using a .pgo file, where we will calculate the
energy levels and transitions 
"""
function PGopher(name::String, pgo_file::String, nproc::Integer=4)
    levels_df, catalog_df = run_pgopher(pgo_file, nproc)
    save("$name.pgopher.jld2", Dict("levels_df" => levels_df, "catalog_df" => catalog_df, "name" => name))
    hash = hash_file("$name.pgopher.jld2")
    combined = merge_pgopher_results(levels_df, catalog_df)
    levels, transitions = dataframes_to_objects(levels_df, combined)
    return PGopher(transitions, levels, name, hash)
end

function PGopher(bson_file::String)
    data_dict = load(bson_file)
    hash = hash_file(bson_file)
    combined = merge_pgopher_results(data_dict["levels_df"], data_dict["catalog_df"]) 
    levels, transitions = dataframes_to_objects(data_dict["levels_df"], combined)
    return PGopher(transitions, levels, data_dict["name"], hash)
end
