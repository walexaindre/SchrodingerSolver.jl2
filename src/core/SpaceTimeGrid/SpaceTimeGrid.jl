include("SpaceTimeGridTypes.jl")

export get_metadata, get_τ, get_measure, get_element, get_diameter
export get_indexed_elements, evaluate_indexed_elements

@inline function check_size(range::AbstractRange{T}, h::T, N::Int) where {T <: AbstractFloat}
    if !(size(range, 1) == N)
    elseif !isapprox(range.step|>T, h)
        throw(ArgumentError("Range step and h are very different... (Possibly miscalculated in SpaceTimeGrid construction) h: $(h), step: $(range.step)"))
    end
end

@inline function SpaceTimeGrid1D(xrange::AbstractRange{T},
    τ::T,
    hx::T,
    Nx::Int) where {T <: AbstractFloat}
    #checking step
    check_size(xrange, hx, Nx)
    SpaceTimeGrid1D(xrange, MetaMesh1D(Nx), τ, hx)
end

@inline function SpaceTimeGrid2D(xrange::AbstractRange{T},
    yrange::AbstractRange{T},
    τ::T,
    hx::T,
    hy::T,
    Nx::Int,
    Ny::Int) where {T <: AbstractFloat}
    #checking step
    check_size(xrange, hx, Nx)
    check_size(yrange, hy, Ny)
    SpaceTimeGrid2D(xrange, yrange, MetaMesh2D(Nx, Ny), τ, hx, hy)
end

@inline function SpaceTimeGrid3D(xrange::AbstractRange{T},
    yrange::AbstractRange{T},
    zrange::AbstractRange{T},
    τ::T,
    hx::T,
    hy::T,
    hz::T,
    Nx::Int,
    Ny::Int,
    Nz::Int) where {T <: AbstractFloat}
    #checking step
    check_size(xrange, hx, Nx)
    check_size(yrange, hy, Ny)
    check_size(zrange, hz, Nz)
    SpaceTimeGrid3D(xrange, yrange, zrange, MetaMesh3D(Nx, Ny, Nz), τ, hx, hy, hz)
end

@inline function CreateGrid(xrange::AbstractRange{T},
    τ::T,
    hx::T,
    Nx::Int) where {T <: AbstractFloat}
    SpaceTimeGrid1D(xrange, τ, hx, Nx)
end
@inline function CreateGrid(xrange::AbstractRange{T},
    yrange::AbstractRange{T},
    τ::T,
    hx::T,
    hy::T,
    Nx::Int,
    Ny::Int) where {T <: AbstractFloat}
    SpaceTimeGrid2D(xrange, yrange, τ, hx, hy, Nx, Ny)
end
@inline function CreateGrid(xrange::AbstractRange{T},
    yrange::AbstractRange{T},
    zrange::AbstractRange{T},
    τ::T,
    hx::T,
    hy::T,
    hz::T,
    Nx::Int,
    Ny::Int,
    Nz::Int) where {T <: AbstractFloat}
    SpaceTimeGrid3D(xrange, yrange, zrange, τ, hx, hy, hz, Nx, Ny, Nz)
end

@inline function get_τ(Grid::SpaceTimeGrid)
    Grid.τ
end

@inline function get_measure(Grid::SpaceTimeGrid1D)
    Grid.hx
end

@inline function get_measure(Grid::SpaceTimeGrid2D)
    Grid.hx * Grid.hy
end

@inline function get_measure(Grid::SpaceTimeGrid3D)
    Grid.hx * Grid.hy * Grid.hz
end

@inline function get_square_root_measure(Grid::SpaceTimeGrid)
    sqrt(get_measure(Grid))
end

@inline function get_element(Grid::SpaceTimeGrid1D, index::Int)
    Grid.Ωx[index]
end

@inline function get_element(Grid::SpaceTimeGrid2D, index_x::Int, index_y::Int)
    Grid.Ωx[index_x], Grid.Ωy[index_y]
end

@inline function get_element(Grid::SpaceTimeGrid3D,
    index_x::Int,
    index_y::Int,
    index_z::Int)
    Grid.Ωx[index_x], Grid.Ωy[index_y], Grid.Ωy[index_z]
end

@inline function get_diameter(Grid::SpaceTimeGrid1D)
    return Grid.hx
end

@inline function get_diameter(Grid::SpaceTimeGrid2D)
    return sqrt(Grid.hx^2 + Grid.hy^2)
end

@inline function get_diameter(Grid::SpaceTimeGrid3D)
    return sqrt(Grid.hx^2 + Grid.hy^2 + Grid.hz^2)
end

@inline function get_metadata(Grid::SpaceTimeGrid)
    Grid.metadata
end

@inline function get_space_steps(Grid::SpaceTimeGrid1D)
    Grid.hx
end

@inline function get_space_steps(Grid::SpaceTimeGrid2D)
    Grid.hx, Grid.hy
end

@inline function get_space_steps(Grid::SpaceTimeGrid3D)
    Grid.hx, Grid.hy, Grid.hz
end

###

@inline function linear_indexing(Grid::SpaceTimeGrid)
    meta = get_metadata(Grid)
    max_size = linear_size(meta)

    linear_index_wrap(index) = linear_indexing(index, meta)

    indexes = collect(1:max_size)

    linear_index_wrap.(indexes)
end

@inline function get_indexed_elements(Grid::SpaceTimeGrid)
    @inline function get_element_wrap(index)
        get_element(Grid, index...)
    end
    indexing = linear_indexing(Grid)
    get_element_wrap.(indexing)
end

@inline function evaluate_indexed_elements(Grid::SpaceTimeGrid, fun::Function)
    indexed_collection = get_indexed_elements(Grid)

    @inline function evaluate_wrap(tuple)
        fun(tuple...)
    end

    evaluate_wrap.(indexed_collection)
end

@inline function evaluate_indexed_elements(Grid::SpaceTimeGrid,
    fun::Function,
    add_argument...)
    indexed_collection = get_indexed_elements(Grid)

    @inline function evaluate_wrap(tuple)
        fun(tuple..., add_argument...)
    end

    evaluate_wrap.(indexed_collection)
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