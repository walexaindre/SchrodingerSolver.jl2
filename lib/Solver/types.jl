export PDESolver

struct MetaKernel
    linearOperatorB::LinearOperator
    linearOperatorC::LinearOperator
    preconditionerB::LinearOperator
end

struct PDESolver{Backend<:AbstractBackend,MetaGrid<:MetaMesh}
    Grid::MetaGrid
    Kernel::Array{MetaKernel}
end