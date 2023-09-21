include("Backends/Backends.jl")
export solve



function solve(::Type{Backend},PDE::SchrodingerPDE,Grid::SpaceTimeGrid,time_order::Int,space_order::Int) where Backend<:AbstractBackend
    
end