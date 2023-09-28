struct SymmetricTimeCompositionMethod{T<:AbstractFloat,C<:AbstractArray{T}}
    order::Int64
    substeps::Int64
    coefficients::C
end