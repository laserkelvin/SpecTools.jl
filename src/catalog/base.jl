

export 
    Catalog,
    Molecule, 
    Q!


abstract type Catalog end

mutable struct Molecule
    name::String
    column_density::AbstractFloat
    temperature::AbstractFloat
    Q::AbstractFloat
    catalog::Catalog

    function Molecule(name::String, catalog::Catalog)
        new(name, 1., 300., 1., catalog)
    end
end

"""
Update the partition function of a molecule object, given a
target temperature.
"""
function Q!(molecule::Molecule, T::AbstractFloat)
    molecule.Q = Q(molecule.levels, T)
end


include("pgopher.jl")
include("spcat.jl")