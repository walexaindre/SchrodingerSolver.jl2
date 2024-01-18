abstract type SolverOptions end

abstract type SolverOptionsGPU <: SolverOptions end

abstract type SolverOptionsCPU <: SolverOptions end

abstract type SolverOptionsCPUParallel <: SolverOptions end

mutable struct SchrodingerSolverOptions{
    ToleranceType <: AbstractFloat,
    Backend <: AbstractBackend,
}
    is_initialized::Bool
    r_tol::ToleranceType
    a_tol::ToleranceType
    fixed_innerloop_steps::Int
    showprogress::Bool
    verbose::Int
    stats::PDESolverStats
    PDE::SchrodingerPDE
    Solver::PDESolver2
    compute_backend::Type{Backend}
    data_type::Type{ToleranceType}
end

SolverOptionsGPUF64 = SchrodingerSolverOptions{Float64, GPUBackend} <: SolverOptionsGPU
SolverOptionsGPUF32 = SchrodingerSolverOptions{Float32, GPUBackend} <: SolverOptionsGPU

SolverOptionsCPUF64 = SchrodingerSolverOptions{Float64, CPUBackend} <: SolverOptionsCPU
SolverOptionsCPUF32 = SchrodingerSolverOptions{Float32, CPUBackend} <: SolverOptionsCPU

SolverOptionsCPUParallelF64 = SchrodingerSolverOptions{Float64, CPUParallelBackend} <:
                              SolverOptionsCPUParallel
SolverOptionsCPUParallelF32 = SchrodingerSolverOptions{Float32, CPUParallelBackend} <:
                              SolverOptionsCPUParallel