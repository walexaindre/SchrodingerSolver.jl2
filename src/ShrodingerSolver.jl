module ShrodingerSolver

using Preferences
using Base.Threads
using Printf
using LinearAlgebra

using CUDA
CUDA.allowscalar(false)

using Krylov
using LinearOperators
using SparseArrays
using GLMakie

using IncompleteLU

export compute_backend,set_compute_backend

include("../lib/index.jl")


end
