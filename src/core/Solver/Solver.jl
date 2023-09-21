include("Backends/Backends.jl")
export solve



function solve(::Type{Backend},PDE::SchrodingerPDE,Grid::SpaceTimeGrid) where Backend<:AbstractBackend
    
end