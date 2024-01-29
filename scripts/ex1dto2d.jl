using SchrodingerSolver
using CUDA
using GLMakie
CUDA.device!(1)
CUDA.allowscalar(false)

boundaries = [0 0; 2*pi 2*pi]
T = 10.0
σ= [1.0;1.0]

function ψ01(x,y,t=0)
    exp(1im*(2*x+y-t))
end


function ψ02(x,y,t=0)
    exp(1im*(2*x+2*y-t))
end

## f1 = x+2/3 y+z
## f2 = 2/3 x+y+2/3 z
## f3 = x+ 2/3 y + z


#  PDE_Meta.N(prevIteration,nextIteration,index)

function N0(prev, current, index::Int)
    p_idx = view(prev, :, index)

    out = (abs2.(p_idx) + abs2.(current))
    if index == 2
        out .*=-3.5
    elseif index == 1
        out.*=-2
    end
    out
end

function FieldF(A)
    A = abs2.(A)
    0.5*(7.0*A[:,2].^2+4.0*A[:,1].^2)
end

Start=[ψ01,ψ02]

PDE = SchrodingerPDEPolynomic(boundaries, σ, N0, Start, T, FieldF)
#@show (2*pi/80)/(8*pi)
#@pprof
#solve(Float64,CPUBackend,PDE,4,τ=(2*pi/400)/(8*pi),Nx=400,Ny=400)
#sol,mesh= solve(Float64,GPUBackend,PDE,2,τ=(2*pi/100)/(8*pi),Nx=200,Ny=200);
#sol,mesh= solve(Float64,CPUBackend,PDE,2,τ=(2*pi/500)/(8*pi),Nx=500,Ny=500);

Opt,Memory = initialize(Float64,GPUBackend,PDE,2,τ=(2*pi/100)/(8*pi),Nx=200,Ny=200);

yval  = Observable(abs2.(Memory.current_state[:,1])|>Array)

surf = surface(Opt.Solver.Grid[:,1],Opt.Solver.Grid[:,2],yval,axis=(type=Axis3, azimuth = pi/4))


display(surf)

for i in 1:5                                                                                                                                                                                                         
    step!(Opt,Memory)                                                                                                                                                                                                               
    Error(Opt,Memory)|>println
    yval[] =   abs2.(Memory.current_state[:,1])|>Array                                                                                               
end