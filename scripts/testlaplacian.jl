using SchrodingerSolver
using LinearAlgebra
using SparseArrays
#using Cthulhu


τ = 0.5
hx=2.0
hy=1.0

xran = 1.0:hx:10.0

yran = 1.0:hy:4.0

Nx = length(xran)
Ny = length(yran)

Grid2D = SpaceTimeGrid2D(xran,yran,τ,hx,hy,Nx,Ny)
M2D = get_metadata(Grid2D)

Grid1D = SpaceTimeGrid1D(xran,τ,hx,Nx)
M1D = get_metadata(Grid1D)
Grid1D2 = SpaceTimeGrid1D(xran,τ,hy,Ny)
M1D2 = get_metadata(Grid1D2)

A = get_sparsematrix_A(Float64,M2D,2)
D = get_sparsematrix_D(Grid2D,2)
display(D|>Array)

A1D = get_sparsematrix_A(Float64,M1D,2)
D1D = get_sparsematrix_D(Grid1D,2)

A1D2 = get_sparsematrix_A(Float64,M1D2,2)
D1D2 = get_sparsematrix_D(Grid1D2,2)

Ix = sparse((1.0)*I, Nx, Nx)
Iy = sparse((1.0)*I, Ny, Ny)

(kron(Iy,D1D)+kron(D1D2,Ix))==D
kronsum()