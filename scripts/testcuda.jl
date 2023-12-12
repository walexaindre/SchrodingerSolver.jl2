using SparseArrays, LinearAlgebra, CUDA, IncompleteLU
using CUDA.CUSPARSE
using CUDA.CUSOLVER

CUDA.allowscalar(false)
using IncompleteLU: ILUFactorization
import LinearAlgebra: ldiv!




A_cpu = sprand(Float32, 1000, 1000, 10 / 1000) + 100I
Precilu = ilu(A_cpu, τ = 3.0)

#LuGPU = Precilu|>CuSparseMatrixCSC


#z = CUDA.zeros(Float32,1000)


function ldiv_ilu0!(P::CuSparseMatrixCSC, x, y, z)
    ldiv!(z, LowerTriangular(P), x)      # Forward substitution with L
    ldiv!(y, UnitUpperTriangular(P), z)  # Backward substitution with U
    return y
    end


#ldiv_ilu0!(Lu)
Precilu

