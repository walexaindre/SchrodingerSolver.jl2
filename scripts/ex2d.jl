using SchrodingerSolver
using PProf
boundaries = [0 0; 2*pi 2*pi]
T = 10.0
σ= [1.0/2,1.0/2,1.0/2.0]

ω₁ = 31.0/6.0
ω₂ = 19.0/3.0
ω₃ = ω₁

function ψ01(x,y,t=0)
    exp(1im*(x+2*y-ω₁*t))
end


function ψ02(x,y,t=0)
    exp(1im*(2*x+2*y-ω₁*t))
end


function ψ03(x,y,t=0)
    exp(1im*(2*x+y-ω₃*t))
end

## f1 = x+2/3 y+z
## f2 = 2/3 x+y+2/3 z
## f3 = x+ 2/3 y + z


#  PDE_Meta.N(prevIteration,nextIteration,index)


function N0(prev,current,index::Int)     
    p_idx = view(prev,:,index)
    out = 0.5*(abs2.(p_idx)+abs2.(current))
    if index==1
        v2 = abs2.(view(prev,:,2))
        v3 = abs2.(view(prev,:,3))
        out+= (2.0/3.0)*v2+v3
    elseif index==2
        v1 = abs2.(view(prev,:,1))
        v3 = abs2.(view(prev,:,3))
        out+= (2.0/3.0)*(v1.+v3)
    elseif index==3
        v1 = abs2.(view(prev,:,1))
        v2 = abs2.(view(prev,:,2))
        out+= v1.+(2.0/3.0)*v2
    end
    out
end

function FieldF(A)
    0.5*(A[:,1].^2+A[:,2].^2+A[:,3].^2)+(2.0/3.0)*(A[:,1].*A[:,2]+A[:,2].*A[:,3])+A[:,1].*A[:,3]
end

Start=[ψ01,ψ02,ψ03]

PDE = SchrodingerPDEPolynomic(boundaries, σ, N0, Start, T, FieldF)
#@show (2*pi/80)/(8*pi)
#@pprof
solve(Float64,CPUBackend,PDE,2,τ=(2*pi/100)/(8*pi),Nx=100,Ny=100)