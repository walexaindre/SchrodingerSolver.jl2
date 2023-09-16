using SchrodingerSolver
using CUDA
using CUDA.CUSPARSE


using SparseArrays
using IncompleteLU

z = MetaMesh2D(4000,6000)

B = generate_kernel_B(Float64,z,1.0,1.0,1.0,1.0)

C = generate_kernel_C(Float64,z,1.0,1.0,1.0,1.0)

PreB = generate_preconditioner_B(B)

meta = assembly_metakernel(CPUBackend,z,B,C,PreB)
test = ones(ComplexF64,24000000)


function calloc()
    meta.preconditionerB*test
end

@time calloc()