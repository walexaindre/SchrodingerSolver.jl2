include("SpaceTimeGridTypes.jl")

export get_metadata, get_τ, get_measure, get_element, get_diameter

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

function get_element(Grid::SpaceTimeGrid3D, index_x::Int, index_y::Int, index_z::Int)
    Grid.Ωx[index_x], Grid.Ωy[index_y], Grid.Ωy[index_z]
end

function get_diameter(Grid::SpaceTimeGrid1D)
    return Grid.hx
end

function get_diameter(Grid::SpaceTimeGrid2D)
    return sqrt(Grid.hx^2 + Grid.hy^2)
end

function get_diameter(Grid::SpaceTimeGrid3D)
    return sqrt(Grid.hx^2 + Grid.hy^2 + Grid.hz^2)
end

function get_metadata(Grid::SpaceTimeGrid)
    Grid.metadata
end

#Base functionality

@inline Base.eltype(Grid::SpaceTimeGrid1D) = eltype(Grid.Ωx)
@inline Base.eltype(Grid::SpaceTimeGrid2D) = (eltype(Grid.Ωx), eltype(Grid.Ωx))
@inline function Base.eltype(Grid::SpaceTimeGrid3D)
    (eltype(Grid.Ωx), eltype(Grid.Ωx), eltype(Grid.Ωx))
end

@inline Base.length(Grid::SpaceTimeGrid) = length(get_metadata(Grid))

@inline Base.size(Grid::SpaceTimeGrid) = size(get_metadata(Grid))

@inline Base.ndims(Grid::SpaceTimeGrid) = ndims(get_metadata(Grid))

@inline Base.firstindex(Grid::SpaceTimeGrid) = firstindex(get_metadata(Grid))

@inline Base.lastindex(Grid::SpaceTimeGrid) = lastindex(get_metadata(Grid))

@inline function Base.iterate(Grid::SpaceTimeGrid1D, state::Int = 1)
    if state > length(Grid)
        return nothing
    else
        x = linear_indexing(state, get_metadata(Grid))
        return (Grid.Ωx[x], state + 1)
    end
end

@inline function Base.iterate(Grid::SpaceTimeGrid2D, state::Int = 1)
    if state > length(Grid)
        return nothing
    else
        x, y = linear_indexing(state, get_metadata(Grid))
        return ((Grid.Ωx[x], Grid.Ωy[y]), state + 1)
    end
end

@inline function Base.iterate(Grid::SpaceTimeGrid3D, state::Int = 1)
    if state > length(Grid)
        return nothing
    else
        x, y, z = linear_indexing(state, get_metadata(Grid))
        return ((Grid.Ωx[x], Grid.Ωy[y], Grid.Ωz[z]), state + 1)
    end
end

@inline function Base.getindex(Grid::SpaceTimeGrid1D, index::Int)
    @boundscheck begin
        if !(1 <= index <= length(Grid.metadata))
            throw(BoundsError(Grid, index))
        end
    end
    Grid.Ωx[linear_indexing(index, Grid.metadata)]
end

@inline function Base.getindex(Grid::SpaceTimeGrid2D, index::Int)
    @boundscheck begin
        if !(1 <= index <= length(Grid.metadata))
            throw(BoundsError(Grid, index))
        end
    end
    x, y = linear_indexing(index, Grid.metadata)
    Grid.Ωx[x], Grid.Ωy[y]
end

@inline function Base.getindex(Grid::SpaceTimeGrid3D, index::Int)
    @boundscheck begin
        if !(1 <= index <= length(Grid.metadata))
            throw(BoundsError(Grid, index))
        end
    end
    x, y, z = linear_indexing(index, Grid.metadata)
    Grid.Ωx[x], Grid.Ωy[y], Grid.Ωz[z]
end