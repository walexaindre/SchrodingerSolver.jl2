using SparseArrays, LinearAlgebra, CUDA, IncompleteLU, Krylov,LinearOperators
using CUDA.CUSPARSE
using IncompleteLU: ILUFactorization
import LinearAlgebra: ldiv!
CUDA.allowscalar(false)

if CUDA.functional()
      # CPU Arrays
  A_cpu = sprand(2000, 2000, 0.003)+300I
  b_cpu = rand(2000)


  b_gpu = CuVector(b_cpu)

  # Transfer the linear system from the CPU to the GPU
  A_gpu = CuSparseMatrixCSR(A_cpu)  # A_gpu = CuSparseMatrixCSC(A_cpu)
  b_gpu = CuVector(b_cpu)

  # ILU(0) decomposition LU ≈ A for CuSparseMatrixCSC or CuSparseMatrixCSR matrices
  P = ilu02(A_gpu)

  # Additional vector required for solving triangular systems
  n = length(b_gpu)
  T = eltype(b_gpu)
  z = CUDA.zeros(T, n)

  # Solve Py = x
  function ldiv_ilu0!(P::CuSparseMatrixCSR, x, y, z)
    ldiv!(z, UnitLowerTriangular(P), x)  # Forward substitution with L
    ldiv!(y, UpperTriangular(P), z)      # Backward substitution with U
    return y
  end

  function ldiv_ilu0!(P::CuSparseMatrixCSC, x, y, z)
    ldiv!(z, LowerTriangular(P), x)      # Forward substitution with L
    ldiv!(y, UnitUpperTriangular(P), z)  # Backward substitution with U
    return y
  end

  x = CUDA.ones(T,n)
  y = CUDA.zeros(T,n)

  ldiv_ilu0!(P, x, y, z)
  P
end



struct CudaILUFactorization{TL, TU}
    L::TL
    U::TU
end

function CudaILUFactorization(f::ILUFactorization)
    L = UnitLowerTriangular(CuSparseMatrixCSR(f.L))
    U = UpperTriangular(CuSparseMatrixCSR(transpose(f.U)))

    CuSparseMatrixCSR()

    return CudaILUFactorization(L, U)
end
