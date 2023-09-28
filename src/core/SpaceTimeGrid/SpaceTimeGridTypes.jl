export SpaceTimeGrid,SpaceTimeGrid1D,SpaceTimeGrid2D,SpaceTimeGrid3D
abstract type SpaceTimeGrid end

struct SpaceTimeGrid1D{DataType<:AbstractFloat} <: SpaceTimeGrid
    Ωx::AbstractRange{DataType}

    metadata::MetaMesh1D

    τ::DataType
    hx::DataType
end


struct SpaceTimeGrid2D{DataType<:AbstractFloat} <: SpaceTimeGrid

    Ωx::AbstractRange{DataType}
    Ωy::AbstractRange{DataType}

    metadata::MetaMesh2D

    τ::DataType
    hx::DataType
    hy::DataType
end


struct SpaceTimeGrid3D{DataType<:AbstractFloat} <: SpaceTimeGrid

    Ωx::AbstractRange{DataType}
    Ωy::AbstractRange{DataType}
    Ωz::AbstractRange{DataType}

    metadata::MetaMesh3D

    τ::DataType
    hx::DataType
    hy::DataType
    hz::DataType
end