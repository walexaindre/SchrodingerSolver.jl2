include("types.jl")

function get_τ(Grid::SpaceTimeGrid)
    Grid.τ
end

function get_measure(Grid::SpaceTimeGrid1D)
    Grid.hx
end

function get_measure(Grid::SpaceTimeGrid2D)
    Grid.hx * Grid.hy
end

function get_measure(Grid::SpaceTimeGrid3D)
    Grid.hx * Grid.hy * Grid.hz
end

function get_element(Grid::SpaceTimeGrid1D, index::Int)
    Grid.Ωx[index]
end

function get_element(Grid::SpaceTimeGrid2D, index_x::Int, index_y::Int)
    Grid.Ωx[index_x], Grid.Ωy[index_y]
end

function get_element(Grid::SpaceTimeGrid3D, index_x::Int, index_y::Int,index_z::Int)
    Grid.Ωx[index_x], Grid.Ωy[index_y],Grid.Ωy[index_z]
end

function get_diameter(Grid::SpaceTimeGrid1D)
    return Grid.hx
end

function get_diameter(Grid::SpaceTimeGrid2D)
    return sqrt(Grid.hx^2+Grid.hy^2)   
end

function get_diameter(Grid::SpaceTimeGrid3D)
    return sqrt(Grid.hx^2+Grid.hy^2+Grid.hz^2)    
end

