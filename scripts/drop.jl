using SchrodingerSolver
using CUDA
using SparseArrays
using Base.Threads
using ProgressMeter
using Krylov

CUDA.device!(1)
CUDA.allowscalar(false)

using PProf
boundaries = [0 0; 2*pi 2*pi]
T = 10.0
σ = [1.0 / 2, 1.0 / 2, 1.0 / 2.0]

ω₁ = 31.0 / 6.0
ω₂ = 19.0 / 3.0
ω₃ = ω₁

function ψ01(x, y, t = 0)
    exp(1im * (x + 2 * y - ω₁ * t))
end

function ψ02(x, y, t = 0)
    exp(1im * (2 * x + 2 * y - ω₁ * t))
end

function ψ03(x, y, t = 0)
    exp(1im * (2 * x + y - ω₃ * t))
end

## f1 = x+2/3 y+z
## f2 = 2/3 x+y+2/3 z
## f3 = x+ 2/3 y + z

#  PDE_Meta.N(prevIteration,nextIteration,index)

function N0(prev, current, index::Int)
    p_idx = view(prev, :, index)

    out = -0.5 * (abs2.(p_idx) - abs2.(current))
    if index == 1
        v2 = abs2.(view(prev, :, 2))
        v3 = abs2.(view(prev, :, 3))
        out -= (2.0 / 3.0) * v2 + v3
    elseif index == 2
        v1 = abs2.(view(prev, :, 1))
        v3 = abs2.(view(prev, :, 3))
        out -= (2.0 / 3.0) * (v1 .+ v3)
    elseif index == 3
        v1 = abs2.(view(prev, :, 1))
        v2 = abs2.(view(prev, :, 2))
        out -= v1 .+ (2.0 / 3.0) * v2
    end
    out
end

function FieldF(A)
    0.5 * (A[:, 1] .^ 2 + A[:, 2] .^ 2 + A[:, 3] .^ 2) +
    (2.0 / 3.0) * (A[:, 1] .* A[:, 2] + A[:, 2] .* A[:, 3]) + A[:, 1] .* A[:, 3]
end

Start = [ψ01, ψ02, ψ03]

PDE = SchrodingerPDEPolynomic(boundaries, σ, N0, Start, T, FieldF)
#@show (2*pi/80)/(8*pi)
#@pprof
#solve(Float64,CPUBackend,PDE,4,τ=(2*pi/400)/(8*pi),Nx=400,Ny=400)
B = solve(Float64, GPUBackend, PDE, 4, τ = (2 * pi / 200) / (8 * pi), Nx = 300, Ny = 200)

function drop(M::SparseMatrixCSC, τ = 0.0001::AbstractFloat)
    #S[I(k),J(k)]=V(k)
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{ComplexF64}()

    lockIJV = ReentrantLock()

    sizehint!(I, ceil(Int64, nnz(M) * 1.3))
    sizehint!(J, ceil(Int64, nnz(M) * 1.3))
    sizehint!(V, ceil(Int64, nnz(M) * 1.3))

    rows, cols = size(M)
    p= Progress(cols)
    #@threads 
    for col_idx in 1:cols
        eₖ = zeros(ComplexF64, rows)
        eₖ[col_idx:col_idx] .= 1.0

        col, _ = gmres(M, eₖ)

        bitindex = abs2.(col) .> τ
        int_index = findall(bitindex)
        lock(lockIJV) do
            append!(I, int_index)
            append!(J, fill(col_idx, length(int_index)))
            append!(V, col[bitindex])
        end
        next!(p)

        if col_idx==5
            break
        end 
    end

    sparse(I, J, V,rows,cols)
end

function dropgpu(M::CUSPARSE.CuSparseMatrixCSR, τ = 0.0001::AbstractFloat)

    

    #S[I(k),J(k)]=V(k)
    I = CUDA.zeros(Int32,ceil(Int32, nnz(M) * 1.3))
    J = CUDA.zeros(Int32,ceil(Int32, nnz(M) * 1.3))
    V = CUDA.zeros(eltype(M),ceil(Int32, nnz(M) * 1.3))

    lockIJV = ReentrantLock()
    current_idx = 0

    rows, cols = size(M)
    
    p= Progress(cols)
    for col_idx in 1:cols
        eₖ = CUDA.zeros(eltype(M), rows)
        eₖ[col_idx:col_idx] .= 1.0
        col, _ = gmres(M, eₖ)

        bitindex = abs2.(col) .> τ

        int_index = findall(bitindex)
        
        itlen = length(int_index)

        lock(lockIJV) do
            I[current_idx+1:current_idx+itlen].=int_index
            J[current_idx+1:current_idx+itlen].=col_idx
            V[current_idx+1:current_idx+itlen].=col[bitindex]
            current_idx+=itlen
        end
        next!(p)
    end
    I = I[begin:current_idx],J[begin:current_idx],V[begin:current_idx]
end