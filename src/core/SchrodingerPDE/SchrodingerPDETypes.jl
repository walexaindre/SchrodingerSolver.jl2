export SchrodingerPDE,SchrodingerPDENonPolynomic,SchrodingerPDEPolynomic
abstract type SchrodingerPDE end

"Generic handler for non polynomic potentials"
struct SchrodingerPDENonPolynomic{DataType<:AbstractFloat,Fun1<:Function,Fun2<:Function} <: SchrodingerPDE
    boundaries::AbstractArray{DataType}
    σ::AbstractArray{DataType,2}
    f::Array{Fun1,1} #Component function
    ψ₀::Array{Fun2,1} #Initial condition
    T::DataType #Time end
    F::Function #field
end



"Optimized structure for polynomical potentials"
struct SchrodingerPDEPolynomic{DataType<:AbstractFloat,Fun1<:Function,Fun2<:Function} <: SchrodingerPDE
    boundaries::AbstractArray{DataType}
    σ::AbstractArray{DataType,2}
    N::Fun1 #Component function optimization
    ψ₀::Array{Fun2,1} #Initial condition
    T::DataType #Time end
    F::Function #field
end