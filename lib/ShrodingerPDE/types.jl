abstract type ShrodingerPDEBase end

"Generic handler for non polynomic potentials"
struct ShrodingerPDENonPolynomic{DataType<:AbstractFloat,Fun1<:Function,Fun2<:Function} <: ShrodingerPDEBase
    boundaries::AbstractArray{DataType}
    σ::AbstractArray{DataType,1}
    f::Array{Fun1,1} #Component function
    ψ₀::Array{Fun2,1} #Initial condition
    T::DataType #Time end
    F::Function #field
end



"Optimized structure for polynomical potentials"
struct ShrodingerPDEPolynomic{DataType<:AbstractFloat,Fun1<:Function,Fun2<:Function} <: ShrodingerPDEBase
    boundaries::AbstractArray{DataType}
    σ::AbstractArray{DataType,1}
    N::Fun1 #Component function optimization
    ψ₀::Array{Fun2,1} #Initial condition
    T::DataType #Time end
    F::Function #field
end