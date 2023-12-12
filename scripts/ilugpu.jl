using SchrodingerSolver
using CUDA
using LinearAlgebra
using SparseArrays
using LinearOperators
using IncompleteLU
using Krylov
using CUDA.CUSPARSE
CUDA.device!(1)
CUDA.allowscalar(false)


Mesh2D = MetaMesh2D(100,100)

Grid = SpaceTimeGrid2D(1.0:400.0,1.0:500.0,Mesh2D,0.0025,0.06346651825433926,0.06346651825433926)
D= 4im*get_sparsematrix_A(Float64,Mesh2D,2) - 0.0025*get_sparsematrix_D(Grid,2)


cpu = ilu(D)

Upper = cpu.U |> CuSparseMatrixCSR
vec = CUDA.ones(ComplexF64,size(D,1))



function backward_substitutionx!(U, y::AbstractVector)

    @inbounds for col = U.dims[1] : -1 : 1

        # Substitutions
        for idx = U.rowPtr[col + 1] - 1 : -1 : U.rowPtr[col] + 1
            y[col] -= U.nzVal[idx] * y[U.colVal[idx]]
        end
        
        # Final answer for y[col]
        y[col] /= U.nzVal[U.rowPtr[col]]
    end

    y
end

backward_substitutionx!(Upper,vec)