export SpaceTimeGrid,SpaceTimeGrid1D,SpaceTimeGrid2D,SpaceTimeGrid3D
abstract type SpaceTimeGrid end

struct SpaceTimeGrid1D{DataType<:AbstractFloat,MetadataMesh<:MetaMesh} <: SpaceTimeGrid
    Ωx::AbstractArray{DataType,1}

    metadata::MetadataMesh

    τ::DataType
    hx::DataType
end


struct SpaceTimeGrid2D{DataType<:AbstractFloat,MetadataMesh<:MetaMesh} <: SpaceTimeGrid

    Ωx::AbstractArray{DataType,1}
    Ωy::AbstractArray{DataType,1}

    metadata::MetadataMesh

    τ::DataType
    hx::DataType
    hy::DataType
end


struct SpaceTimeGrid3D{DataType<:AbstractFloat,MetadataMesh<:MetaMesh} <: SpaceTimeGrid

    Ωx::AbstractArray{DataType,1}
    Ωy::AbstractArray{DataType,1}
    Ωz::AbstractArray{DataType,1}

    metadata::MetadataMesh

    τ::DataType
    hx::DataType
    hy::DataType
    hz::DataType
end