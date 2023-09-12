struct MetaKernel
    linearOperatorB::LinearOperator
    linearOperatorC::LinearOperator
    preconditionerB::LinearOperator
    preconditionerC::LinearOperator
end

struct PDESolver{Backend<:AbstractBackend,MetaGrid<:MetaMesh}
    Grid::MetaGrid
    Kernel::Array{MetaKernel}
end