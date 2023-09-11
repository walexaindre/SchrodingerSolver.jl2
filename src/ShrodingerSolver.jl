module ShrodingerSolver

using Preferences
using Base.Threads
using Printf
using LinearAlgebra

using CUDA

using Krylov
using LinearOperators
using SparseArrays
using GLMakie

export compute_backend

include("../lib/index.jl")


end
