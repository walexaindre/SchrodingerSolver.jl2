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



struct LUperso
	L
	Ut	# transpose of U in LU decomposition
end

Precilu=lu(D)
cpu = ilu(D)
ilucpu = cpu.L+cpu.U'
iluhamgpu = ilucpu|>CuSparseMatrixCSR


Mgpu = CUDA.CUSPARSE.CuSparseMatrixCSR(D)
P = CUDA.CUSPARSE.ilu02(Mgpu)

#throw("Error")
#=Precilu_gpu = LUperso(LowerTriangular(CUDA.CUSPARSE.CuSparseMatrixCSR(I+Precilu.L)), UpperTriangular(CUDA.CUSPARSE.CuSparseMatrixCSR(sparse(Precilu.U'))));

function ldiv_ilu0!(LUp::LUperso, x, y, z)
    ldiv!(z, LUp.L, x)      # Forward substitution with L
    ldiv!(y, LUp.Ut, z)  # Backward substitution with U
    return y
    end
=#
    function ldiv_ilu0!(P::CuSparseMatrixCSR, x, y, z)
        ldiv!(z, UnitLowerTriangular(P), x)      # Forward substitution with L
        ldiv!(y, UpperTriangular(P), z)  # Backward substitution with U
        return y
      end

      function ldiv_ilu0!(P::SparseMatrixCSC, x, y, z)
        ldiv!(z, UnitLowerTriangular(P), x)      # Forward substitution with L
        ldiv!(y, UpperTriangular(P), z)  # Backward substitution with U
        return y
      end

      symmetric = hermitian = false
      n=size(Mgpu,1)
      z = CUDA.zeros(ComplexF64,size(D,1))
      opM = LinearOperator(ComplexF64, n, n, symmetric, hermitian, (y, x) -> ldiv_ilu0!(iluhamgpu, x, y, z))
      zcpu = zeros(ComplexF64,size(D,1))
      opMcpu = LinearOperator(ComplexF64, n, n, symmetric, hermitian, (y, x) -> ldiv_ilu0!(ilucpu, x, y, zcpu))


T0 =CUDA.ones(ComplexF64,size(D,1))
T0cpu =rand(ComplexF64,size(D,1))
T1 =CUDA.zeros(ComplexF64,size(D,1))
T1cpu = ones(ComplexF64,size(D,1))

T2 =CUDA.zeros(ComplexF64,size(D,1))
CUDA.@time r = gmres(Mgpu,T0,verbose=1)
@time rr=gmres(D,T0cpu,verbose=1,M=cpu,ldiv=true)
@time Precilu\T0cpu

D*(cpu\T0cpu)
D*(opMcpu*T0cpu)

IncompleteLU.forward_substitution!(zcpu,cpu,T0cpu)
IncompleteLU.backward_substitution!(T1cpu,cpu,zcpu)

D*T1cpu

zcpu

ldiv!(zcpu,UnitLowerTriangular(ilucpu),T0cpu)
ldiv!(T1cpu,UpperTriangular(ilucpu),zcpu)



ldiv!(zcpu,UnitLowerTriangular(cpu.L),T0cpu)

ldiv!(T1cpu,UpperTriangular(ilucpu),zcpu)





opMcpu*T0cpu
cpu\T0cpu

T0cpu

ldiv_ilu0!(ilucpu, T0cpu, T1cpu, zcpu)