
struct Pickett <: BaseCatalog
    transitions::Vector{<:Transition}
    name::String
    hash::String
    
    function Pickett(name::String, file::String)
        hash = bytes2hex(sha2_256(read(file)))
      end
end
