export SolverMethod,PaulMethod,ADIMethod,PDESolution
abstract type SolverMethod end

abstract type PaulMethod <: SolverMethod end

abstract type ADIMethod <: SolverMethod end

OperatorOrMatrix = Union{LinearOperator,AbstractArray,UniformScaling}

struct MetaKernel{T<:OperatorOrMatrix}
    linearOperatorB::T
    linearOperatorC::T
    preconditionerB
end

struct PDESolver{Backend<:AbstractBackend,MetaGrid<:MetaMesh,Operator<:OperatorOrMatrix}
    Grid::MetaGrid
    Kernel::Array{MetaKernel}
    opA::Operator
end

struct PDESolution 
    time_steps::Int #Number of time_steps

end