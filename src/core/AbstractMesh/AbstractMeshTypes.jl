export MetaMesh,MetaMesh1D,MetaMesh2D,MetaMesh3D

abstract type MetaMesh end

struct MetaMesh1D <: MetaMesh
    N::Int #points
end

struct MetaMesh2D <: MetaMesh
    M::Int #rows
    N::Int #cols
    MN::Int #M×N

    MetaMesh2D(M::Int, N::Int) = new(M, N, M * N)
end

struct MetaMesh3D <: MetaMesh
    M::Int #rows
    N::Int #cols
    L::Int #depth
    MN::Int #M×N

    MetaMesh3D(M::Int, N::Int, L::Int) = new(M, N, L, M * N)
end