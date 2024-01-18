#BLAS level optimizations enabled with multithread enabled for BLAS subroutines!
include("CPU/CPU.jl")

#CPU with care of parallelizable things
include("CPUParallel/Preconditioner.jl")
include("CPUParallel/CPUParallel.jl")

#CUDA
include("CUDA/MemoryAlloc.jl")
include("CUDA/Preconditioner.jl")
include("CUDA/CUDA.jl")