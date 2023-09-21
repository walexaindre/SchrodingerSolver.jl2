include("SpaceDiscretizationTypes.jl")
export get_space_discretization_coefficients
export get_sparsematrix_A

function get_space_discretization_coefficients(::Type{T},
    order::Int)::SecondDerivativeCoefficients{T} where {T <: AbstractFloat}
    if order == 2
        α = 0//1
        β = 0//1

        a = 1//1
        b = 0//1
        c = 0//1

    elseif order == 4
        α = 1 // 10
        β = 0 // 1

        a = (4//3)*(1-α)
        b = (1//3)*(-1+10*α)
        c = 0 // 1
    elseif order == 6
        α = 2 // 11
        β = 0 // 1

        a = 12 // 11
        b = 3 // 11
        c = 0 // 1
    elseif order == 8
        α = 344 // 1179
        β = (38 * α - 9) // 214

        a = (696 - 1191 * α) // 428
        b = (2454 * α - 294) // 535
        c = (1179 * α - 344) // 2140

    elseif order==10
        α = 334 // 899
        β = 43 // 1798

        a = 1065 // 1798
        b = 1038 // 899
        c = 79 // 1798
    else 
        throw(BoundsError([2,4,6,8,10],order))
    end

    SecondDerivativeCoefficients{T}(T(a), T(b), T(c), T(α), T(β))
end

function get_sparsematrix_A(::Type{T},Mesh::MetaMesh,order::Int) where T<:AbstractFloat

    space_discretization = get_space_discretization_coefficients(T,order)
    
    count = 1
    value = [T(1)]
    if (space_discretization.α != 0)
        count+=2
        push!(value,space_discretization.α)
        push!(value,space_discretization.α)
    end
    if (space_discretization.β != 0)
        count+=2
        push!(value,space_discretization.β)
        push!(value,space_discretization.β)
    end
    
    space_usage = count*length(Mesh)

    I = zeros(Int64,space_usage) #row idx
    J = zeros(Int64,space_usage) #column idx
    V = zeros(T,space_usage) #value


    for idx in 1:length(Mesh)

        I[count*(index-1)+1:count*index].= idx
        J[count*(index-1)+1:count*index].= 
        V[count*(index-1)+1:count*index].= value
    end


    sparse(I,J,V)    
end