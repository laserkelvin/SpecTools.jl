
module Catalog

using CSV, DataFrames, DelimitedFiles, BSON
using ..SpecTools


export 
    BaseCatalog,
    Molecule, 
    Q!,
    ν

abstract type 
    BaseCatalog 
end

mutable struct Molecule
    name::String
    Ncol::AbstractFloat
    T::AbstractFloat
    Q::AbstractFloat
    catalog::BaseCatalog

    function Molecule(name::String, catalog::BaseCatalog)
        new(name, 1., 300., 1., catalog)
    end
end

frequencies(molecule::Molecule) = [transition.ν for transition in molecule.catalog.transitions]

"""
Update the partition function of a molecule object, given a
target temperature.
"""
function Q!(molecule::Molecule, T::AbstractFloat)
    molecule.Q = Q(molecule.levels, T)
end

export
    read_pgopher_catalog,
    read_pgopher_levels,
    PGopherLevels,
    PGopher,
    Q,
    run_pgopher,
    run_pgopher_catalog,
    run_pgopher_levels,
    read_pgopher_csv,
    merge_pgopher_results

include("pgopher.jl")

# export Pickett

# include("spcat.jl")

end