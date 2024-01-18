"
    Memory storage for fixed structures
"
struct BackendMemoryCUDA{T <: AbstractFloat}
    current_state::CuArray{Complex{T}, 1}
    component_place::CuArray{Complex{T}, 1}
    b_temp::CuArray{Complex{T}, 1}
    b0_temp::CuArray{Complex{T}, 1}
    zcomp::CuArray{Complex{T}, 1}
    solver_memory::GmresSolver{T, Complex{T}, CuArray{Complex{T}, 1}}
end

function initialize_memory(::Type{T},
        element_count::Int,
        memory_size::Int) where {T <: AbstractFloat}
    mem = CUDA.zeros(T, elem_count)

    BackendMemoryCUDA(mem,
        copy(mem),
        copy(mem),
        copy(mem),
        copy(mem),
        copy(mem),
        GmresSolver(element_count, element_count, memory_size, typeof(mem)))
end