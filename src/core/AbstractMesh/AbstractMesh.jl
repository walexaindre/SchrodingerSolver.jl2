include("AbstractMeshTypes.jl")

export linear_indexing, step, step_x, step_y, step_z, linear_size

@inline function linear_indexing(index::Int, metadata::MetaMesh1D)
    index
end

@inline function linear_indexing(index::Int, metadata::MetaMesh2D)
    col, row = divrem(index - 1, metadata.M)
    row + 1, col + 1
end

@inline function linear_indexing(index::Int, metadata::MetaMesh3D)
    z, remaining_position = divrem(index - 1, metadata.MN)
    y, x = divrem(remaining_position, metadata.M)
    x + 1, y + 1, z + 1
end

@inline function linear_indexing(x::Int, y::Int, z::Int, metadata::MetaMesh3D)
    (z - 1) * metadata.MN + (y - 1) * metadata.M + x
end

@inline function linear_indexing(x::Int, y::Int, metadata::MetaMesh2D)
    (y - 1) * metadata.M + x
end

@inline function step(index::Int, lim::Int)
    r1 = (index - 1 < 1) ? lim : index - 1
    r2 = (index + 1 > lim) ? 1 : index + 1
    (r1, r2)
end

@inline function step(index::Int, lim::Int, step_size::Int)
    offset = index - 1
    r1 = mod(offset - step_size, lim) + 1
    r2 = mod(offset + step_size, lim) + 1
    (r1, r2)
end

@inline function step_x(index::Int, metadata::MetaMesh)
    step(index, metadata.M)
end

@inline function step_y(index::Int, metadata::MetaMesh2D)
    step(index, metadata.N)
end

@inline function step_y(index::Int, metadata::MetaMesh3D)
    step(index, metadata.N)
end

@inline function step_z(index::Int, metadata::MetaMesh3D)
    step(index, metadata.L)
end

@inline function step_x(index::Int, step_size::Int, metadata::MetaMesh)
    step(index, metadata.M, step_size)
end

@inline function step_y(index::Int, step_size::Int, metadata::MetaMesh2D)
    step(index, metadata.N, step_size)
end

@inline function step_y(index::Int, step_size::Int, metadata::MetaMesh3D)
    step(index, metadata.N, step_size)
end

@inline function step_z(index::Int, step_size::Int, metadata::MetaMesh3D)
    step(index, metadata.L, step_size)
end

@inline function linear_size(metadata::MetaMesh1D)
    metadata.M
end

@inline function linear_size(metadata::MetaMesh2D)
    metadata.MN
end

@inline function linear_size(metadata::MetaMesh3D)
    metadata.MN * metadata.L
end

#Base methods
@inline Base.length(AbstractMesh::MetaMesh1D) = AbstractMesh.M
@inline Base.length(AbstractMesh::MetaMesh2D) = AbstractMesh.MN
@inline Base.length(AbstractMesh::MetaMesh3D) = AbstractMesh.MN * AbstractMesh.L

@inline function Base.firstindex(Mesh::MetaMesh)
    1
end

@inline function Base.lastindex(Mesh::MetaMesh)
    length(Mesh)
end

@inline Base.eltype(Mesh::MetaMesh1D) = Int64
@inline Base.eltype(Mesh::MetaMesh2D) = (Int64,Int64)
@inline Base.eltype(Mesh::MetaMesh3D) = (Int64,Int64,Int64)
    
function Base.iterate(AbstractMesh::MetaMesh, state::Int = 1)
    state > length(AbstractMesh) ? nothing : (linear_indexing(state,AbstractMesh),state+1)
end

#1D
@inline function Base.getindex(Mesh::MetaMesh1D, index::Int)
    @boundscheck begin
        if !(1 <= index <= Mesh.M)
            throw(BoundsError(Mesh, index))
        end
    end
    index
end

@inline function Base.size(Mesh::MetaMesh1D)
    Mesh.M
end

@inline function Base.ndims(Mesh::MetaMesh1D)
    1
end

#2D
@inline function Base.getindex(Mesh::MetaMesh2D, index::Int)
    @boundscheck begin
        if !(1 <= index <= Mesh.MN)
            throw(BoundsError(Mesh, index))
        end
    end
    linear_indexing(index, Mesh)
end

@inline function Base.getindex(Mesh::MetaMesh2D, i::Int, j::Int)
    @boundscheck begin
        if !((1 <= i <= Mesh.M) && (1 <= j <= Mesh.N))
            throw(BoundsError(Mesh, (i, j)))
        end
    end
    linear_indexing(i, j, Mesh)
end

@inline function Base.size(Mesh::MetaMesh2D)
    (Mesh.M, Mesh.N)
end

@inline function Base.size(Mesh::MetaMesh2D, dim::Int)
    @boundscheck begin
        if !(1 <= dim <= 2)
            throw(BoundsError(Mesh, dim))
        end
    end

    if dim == 1
        return Mesh.M
    end
    Mesh.N
end

@inline function Base.ndims(Mesh::MetaMesh2D)
    2
end

@inline function Base.firstindex(Mesh::MetaMesh2D, d::Int)
    1
end

@inline function Base.lastindex(Mesh::MetaMesh2D, d::Int)
    if d == 1
        return Mesh.M
    elseif d == 2
        return Mesh.N
    else
        throw(BoundsError(Mesh, d))
    end
end

#3D
@inline function Base.getindex(Mesh::MetaMesh3D, index::Int)
    @boundscheck begin
        if !(1 <= index <= linear_size(Mesh))
            throw(BoundsError(Mesh, index))
        end
    end
    linear_indexing(index, Mesh)
end

@inline function Base.getindex(Mesh::MetaMesh3D, i::Int, j::Int, k::Int)
    @boundscheck begin
        if !((1 <= i <= Mesh.M) && (1 <= j <= Mesh.N) && (1 <= k <= Mesh.L))
            throw(BoundsError(Mesh, (i, j, k)))
        end
    end
    linear_indexing(i, j, k, Mesh)
end

@inline function Base.size(Mesh::MetaMesh3D)
    (Mesh.M, Mesh.N, Mesh.L)
end

@inline function Base.size(Mesh::MetaMesh3D, dim::Int)
    @boundscheck begin
        if !(1 <= dim <= 3)
            throw(BoundsError(Mesh, dim))
        end
    end

    if dim == 1
        return Mesh.M
    elseif dim == 2
        return Mesh.N
    end
    Mesh.L
end

@inline function Base.ndims(Mesh::MetaMesh3D)
    3
end

@inline function Base.firstindex(Mesh::MetaMesh3D, d::Int)
    1
end

@inline function Base.lastindex(Mesh::MetaMesh3D, d::Int)
    if d == 1
        return Mesh.M
    elseif d == 2
        return Mesh.N
    elseif d == 3
        return Mesh.L
    else
        throw(BoundsError(Mesh, d))
    end
end

#################################################################