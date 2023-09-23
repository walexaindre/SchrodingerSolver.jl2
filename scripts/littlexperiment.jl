#using SchrodingerSolver
using LinearAlgebra

xsize = 6
ysize = 5
zsize = 7

xmesh = MetaMesh1D(xsize)
ymesh = MetaMesh1D(ysize)
zmesh = MetaMesh1D(zsize)

xymesh = MetaMesh2D(xsize,ysize)
xyzmesh = MetaMesh3D(xsize,ysize,zsize)

xd1 = get_sparsematrix_A(Float64,xmesh,10)
yd1 = get_sparsematrix_A(Float64,ymesh,10)
zd1 = get_sparsematrix_A(Float64,zmesh,10)

xyd2 = get_sparsematrix_A(Float64,xymesh,10)
xyzd3 = get_sparsematrix_A(Float64,xyzmesh,10)

xId = Matrix(I,xsize,xsize)
yId = Matrix(I,ysize,ysize)
zId = Matrix(I,zsize,zsize)

xyd2alg = kron(yId,xd1)+kron(yd1,xId)

xyzd3alg = kron(zId,yId,xd1)+kron(zId,yd1,xId) + kron(zd1,yId,xId)
iszero(xyd2-xyd2alg)
iszero(xyzd3-xyzd3alg)

