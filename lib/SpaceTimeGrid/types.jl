abstract type SpaceTimeGrid end

struct SpaceTimeGrid1D{DataType<:AbstractFloat,MetadataMesh<:MetaMesh} <: SpaceTimeGrid
    Ωx::AbstractArray{DataType,1}

    metadata::MetadataMesh

    const τ::DataType
    const hx::DataType
end


struct SpaceTimeGrid2D{DataType<:AbstractFloat,MetadataMesh<:MetaMesh} <: SpaceTimeGrid

    Ωx::AbstractArray{DataType,1}
    Ωy::AbstractArray{DataType,1}

    metadata::MetadataMesh

    const τ::DataType
    const hx::DataType
    const hy::DataType
end


struct SpaceTimeGrid3D{DataType<:AbstractFloat,MetadataMesh<:MetaMesh} <: SpaceTimeGrid

    Ωx::AbstractArray{DataType,1}
    Ωy::AbstractArray{DataType,1}
    Ωz::AbstractArray{DataType,1}

    metadata::MetadataMesh

    const τ::DataType
    const hx::DataType
    const hy::DataType
    const hz::DataType
end
