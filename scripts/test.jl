using SparseMatricesCSR
using CUDA
using CUDA.CUSPARSE


V = ones(Float64,7)
I = [1,2,3,4,5,6,2]
J = [1,2,3,4,5,6,1]

A = sparsecsr(I,J,V)

#CuSparseMatrixCSR{Tv, Ti}(rowPtr::CuVector{<:Integer}, colVal::CuVector{<:Integer},
#nzVal::CuVector, dims::NTuple{2,<:Integer})
#function CuSparseMatrixCSR{Tv, Ti}(rowPtr::CuVector{<:Integer}, colVal::CuVector{<:Integer},
#    nzVal::CuVector, dims::NTuple{2,<:Integer}) where {Tv, Ti<:Integer}
#new{Tv, Ti}(rowPtr, colVal, nzVal, dims, length(nzVal))
#end

function CuSparseMatrixCSRs(Matrix::SparseMatrixCSR)
    CuSparseMatrixCSR(Matrix.rowptr|>CuArray,Matrix.colval|>CuArray,Matrix.nzval|>CuArray,(Matrix.m,Matrix.n))
end

CuSparseMatrixCSRs(A)
#=
using ProgressMeter
@showprogress "Comp" showspeed=true barlen=30 for j in 1:100
    @showprogress "V" offset=1 showspeed=true barlen=30 for i in 1:10
        sleep(0.25)
    end
end

=#