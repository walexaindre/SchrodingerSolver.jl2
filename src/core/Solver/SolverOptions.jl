abstract type SolverOptions end

abstract type SolverOptionsGPU <: SolverOptions end

abstract type SolverOptionsCPU <: SolverOptions end

abstract type SolverOptionsCPUParallel <: SolverOptions end

mutable struct SchrodingerSolverOptions{
    FloatType <: AbstractFloat,
    Backend <: AbstractBackend}
    is_initialized::Bool
    r_tol::FloatType
    a_tol::FloatType
    fixed_innerloop_steps::Int
    max_iterations::Int
    showprogress::Bool
    verbose::Int
    stats::PDESolverStats
    PDE::SchrodingerPDE
    Solver::PDESolver3
    start_power::Vector{FloatType}
    start_energy::FloatType
    time_collection::Vector{FloatType}
    compute_backend::Type{Backend}
    data_type::Type{FloatType}
end