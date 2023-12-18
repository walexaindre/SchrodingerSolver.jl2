include("TimeDiscretizationTypes.jl")
####

@inline get_coefficient(CompositionMethod::SymmetricTimeCompositionMethod,index::Int) = CompositionMethod.coefficients[mod(index-1,length(CompositionMethod.coefficients))+1]

@inline function check_coefficients(CompositionMethod::SymmetricTimeCompositionMethod)
    sum_coefficients = 2*sum(CompositionMethod.coefficients)-CompositionMethod.coefficients[end]
    if !isapprox(sum_coefficients,1) 
        throw(DomainError(sum_coefficients,"Must be nearest to one the sum of substeps if you want a valid SymmetricTimeCompositionMethod"))
    end
end

function ConstructSymmetricTimeCompositionMethod(order::Int,substeps::Int,coefficients::AbstractArray{T}) where {T<:AbstractFloat}
    R = SymmetricTimeCompositionMethod(order,substeps,coefficients)
    check_coefficients(R)
    R
end

@inline Base.length(CompositionMethod::SymmetricTimeCompositionMethod) = CompositionMethod.substeps
@inline Base.firstindex(CompositionMethod::SymmetricTimeCompositionMethod) = 1
@inline Base.lastindex(CompositionMethod::SymmetricTimeCompositionMethod) = length(CompositionMethod)

@inline function Base.getindex(CompositionMethod::SymmetricTimeCompositionMethod, index::Int)
    @boundscheck begin
        if !(1 <= index <= length(CompositionMethod))
            throw(BoundsError(1:length(CompositionMethod), index))
        end
    end
    get_coefficient(CompositionMethod,index)
end

@inline Base.eltype(::Type{SymmetricTimeCompositionMethod{T,V}}) where {T,V} = T

function time_discretization_ord_2(::Type{T},
    substeps::Int64,
    index::Int64) where {T <: AbstractFloat}
    return ConstructSymmetricTimeCompositionMethod(2, 1, convert(Array{T}, [1.0]))
end

function get_coefficients(CompositionMethod::SymmetricTimeCompositionMethod)
    CompositionMethod.coefficients
end

function time_discretization_ord_4(::Type{T},
    substeps::Int64,
    index::Int64) where {T <: AbstractFloat}
    if substeps == 5 && index == 1
        γ₁ = (3.0 + sqrt(3.0)) / 6
        γ₂ = (3.0 - sqrt(3.0)) / 6
        γ₃ = -1.0
        return ConstructSymmetricTimeCompositionMethod(4, 5, convert(Array{T}, [γ₁, γ₂, γ₃]))
    elseif substeps == 5 && index == 2
        γ₁ = 0.28
        γ₂ = 0.62546642846767004501
        γ₃ = 1 − 2(γ₁ + γ₂)
        return ConstructSymmetricTimeCompositionMethod(4, 5, convert(Array{T}, [γ₁, γ₂, γ₃]))
    elseif substeps == 7 && index == 1
        γ₁ = 1 / (6 - cbrt(6))
        γ₂ = γ₁
        γ₃ = γ₁
        γ₄ = 1 - 6 * γ₁
        return ConstructSymmetricTimeCompositionMethod(4, 7, convert(Array{T}, [γ₁, γ₂, γ₃, γ₄]))
    else
        throw(BoundsError([(5, 1), (5, 2), (7, 1)], (substeps, index)))
    end
end

function time_discretization_ord_6(::Type{T}, substeps::Int64,
    index::Int64) where {T <: AbstractFloat}
    if substeps == 7 && index == 1
        γ₁ = 0.7845136104775572638
        γ₂ = 0.2355732133593581336
        γ₃ = -1.1776799841788710069
        γ₄ = 1 - 2 * (γ₁ + γ₂ + γ₃)
        return ConstructSymmetricTimeCompositionMethod(6, 7, convert(Array{T}, [γ₁, γ₂, γ₃, γ₄]))
    elseif substeps == 9 && index == 1
        γ₁ = 0.186
        γ₂ = 0.5554970237124783991
        γ₃ = 0.1294669489134753580
        γ₄ = -0.8432656233877346085
        γ₅ = 1 - 2 * (γ₁ + γ₂ + γ₃ + γ₄)
        return ConstructSymmetricTimeCompositionMethod(6, 9, convert(Array{T}, [γ₁, γ₂, γ₃, γ₄, γ₅]))
    elseif substeps == 9 && index == 2
        γ₁ = 0.3921614440073141392
        γ₂ = 0.3325991367893594381
        γ₃ = -0.7062461725576393598
        γ₄ = 0.0822135962935508002
        γ₅ = 0.7985439909348299631
        return ConstructSymmetricTimeCompositionMethod(6, 9, convert(Array{T}, [γ₁, γ₂, γ₃, γ₄, γ₅]))
    else
        throw(BoundsError([(7, 1), (9, 1), (9, 2)], (substeps, index)))
    end
end

function time_discretization_ord_8(::Type{T}, substeps::Int64,
    index::Int64) where {T <: AbstractFloat}
    if substeps == 15 && index == 1
        γ₁ = 0.74167036435061295345
        γ₂ = -0.40910082580003159400
        γ₃ = 0.19075471029623837995
        γ₄ = -0.57386247111608226666
        γ₅ = 0.29906418130365592384
        γ₆ = 0.33462491824529818378
        γ₇ = 0.31529309239676659663
        γ₈ = -0.7968879393529163540
        return ConstructSymmetricTimeCompositionMethod(8,
            15,
            convert(Array{T}, [γ₁, γ₂, γ₃, γ₄, γ₅, γ₆, γ₇, γ₈]))
    else
        throw(BoundsError([(15, 1)], (substeps, index)))
    end
end

function get_time_discretization(::Type{T},
    order::Int,
    substeps::Int,
    index::Int) where {T <: AbstractFloat}
    if order == 2
        return time_discretization_ord_2(T, substeps, index)
    elseif order == 4
        return time_discretization_ord_4(T, substeps, index)
    elseif order == 6
        return time_discretization_ord_6(T, substeps, index)
    elseif order == 8
        return time_discretization_ord_8(T, substeps, index)
    else
        throw(BoundsError([2, 4, 6, 8], order))
    end
end

##########
# Base Method Overload
##########

@inline function Base.iterate(TimeCompositionMethod::SymmetricTimeCompositionMethod, state::Int = 1)
    if state > length(TimeCompositionMethod)
        return nothing
    else
        return (TimeCompositionMethod[state], state + 1)
    end
end