include("types.jl")

function get_kernel_B(Solver::PDESolver,index::Int)
    Solver.Kernel[index].linearOperatorB
end

function get_kernel_C(Solver::PDESolver,index::Int)
    Solver.Kernel[index].linearOperatorC
end