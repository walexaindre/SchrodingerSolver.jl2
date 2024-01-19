abstract type SolverOptions end

abstract type SolverOptionsGPU <: SolverOptions end

abstract type SolverOptionsCPU <: SolverOptions end

abstract type SolverOptionsCPUParallel <: SolverOptions end

mutable struct SchrodingerSolverOptions{
    ToleranceType <: AbstractFloat,
    Backend <: AbstractBackend}
    is_initialized::Bool
    r_tol::ToleranceType
    a_tol::ToleranceType
    fixed_innerloop_steps::Int
    max_iterations::Int
    showprogress::Bool
    verbose::Int
    stats::PDESolverStats
    PDE::SchrodingerPDE
    Solver::PDESolver3
    time_collection::Vector{ToleranceType}
    compute_backend::Type{Backend}
    data_type::Type{ToleranceType}
end