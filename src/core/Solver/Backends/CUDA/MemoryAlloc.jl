"
    Memory storage for fixed structures
"
struct BackendMemoryCUDA{T <: AbstractFloat}
    current_state::CuArray{Complex{T}, 2}
    component_place::CuArray{Complex{T}, 1}
    b_temp::CuArray{Complex{T}, 1}
    b0_temp::CuArray{Complex{T}, 1}
    zcomp::CuArray{Complex{T}, 1}
    solver_memory::GmresSolver{T, Complex{T}}
end

function initialize_memory(::Type{T},ncomponents::Int,
        element_count::Int,
        memory_size::Int) where {T <: AbstractFloat}
    
    
    state = CUDA.zeros(Complex{T}, element_count, ncomponents)
    mem = CUDA.zeros(Complex{T}, element_count)


    BackendMemoryCUDA{T}(
        state,
        mem,
        copy(mem),
        copy(mem),
        copy(mem),
        GmresSolver(element_count, element_count, memory_size, typeof(mem)))
end