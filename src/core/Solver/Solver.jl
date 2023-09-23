include("SolverTypes.jl")
include("Backends/Backends.jl")
export solve

function get_time_steps(PDE::SchrodingerPDE, Grid::SpaceTimeGrid)
    PDE.T / Grid.τ
end




function solve(::Type{Backend},PDE::SchrodingerPDE,Grid::SpaceTimeGrid,time_order::Int,space_order::Int) where Backend<:AbstractBackend
    
end