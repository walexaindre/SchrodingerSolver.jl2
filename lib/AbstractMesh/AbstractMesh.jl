include("types.jl")

function linear_indexing(index::Int, metadata::MetaMesh1D)
    index
end

function linear_indexing(index::Int, metadata::MetaMesh2D)
    col, row = divrem(index - 1, metadata.M)
    row + 1, col + 1
end

function linear_indexing(index::Int, metadata::MetaMesh3D)
    z, remaining_position = divrem(index - 1, metadata.MN)
    y, x = divrem(remaining_position, metadata.M)
    x + 1, y + 1, z + 1
end

function linear_indexing(x::Int, y::Int, z::Int, metadata::MetaMesh3D)
    (z - 1) * metadata.MN + (y - 1) * metadata.M + x
end

function linear_indexing(x::Int, y::Int, metadata::MetaMesh2D)
    (y - 1) * metadata.M + x
end

function step(index::Int, lim::Int)
    r1 = max(index - 1, 1)
    r2 = min(index + 1, lim)
    r1, r2
end

function step_x(index::Int, metadata::MetaMesh)
    step(index, metadata.N)
end

function step_y(index::Int, metadata::MetaMesh2D)
    step(index, metadata.M)
end

function step_y(index::Int, metadata::MetaMesh3D)
    step(index, metadata.M)
end

function step_z(index::Int, metadata::MetaMesh3D)
    step(index, metadata.L)
end

function linear_size(metadata::MetaMesh1D)
    metadata.N
end

function linear_size(metadata::MetaMesh2D)
    metadata.MN
end

function linear_size(metadata::MetaMesh3D)
    metadata.MN * metadata.L
end