using SchrodingerSolver

function ψ01(x,y)
    k1 = 1.0+1.0im
    k1sqr = k1*k1
    l1=k1
    l1sqr = k1sqr
    l1star = 1.0-1.0im
    k1star = l1star

    eta = 0.5*(2.0/(4.0+4.0))

    exp_eta = 0.5*log(eta)

    sqrteta = 1.0/sqrt(eta)
    
    theta = k1*x+l1*y

    0.5*sqrteta*sech(real(theta)+exp_eta)*exp(1.0im*imag(theta))    
end

boundaries = [-20.0 -20.0; 20.0 20.0]

σ = [1.0 1.0]


Start = [ψ01, ψ01]

T = 40.0

α = 1.0

hx=8.0/100
hy=8.0/100

τ = 0.005
function Nop(prevMove,currMove,index::Int)
    #0.5*sum(abs2,currMove,dims=2)+c*view(prevMove,:,3-index)    
    0.5*(abs2.(view(currMove,:,index))+abs2.(view(prevMove,:,index)))+α*abs2.(view(prevMove,:,3-index))   
end

function FieldF(A)
    c₁ = (α |> real(eltype(A)))
    d₁ = (0.5 |> real(eltype(A)))

    B = abs2.(A)
    B₁ = @view B[:, 1]
    B₂ = @view B[:, 2]
    d₁ * (B₁ .^ 2 .+ B₂ .^ 2) .+ c₁ * B₁ .* B₂
end


function ψ1(B)
    A = abs2.(B)
    @views A[:, 1] .+ (α |> eltype(A[1, :])) .* A[:, 2]
end


function ψ2(B)
    A = abs2.(B)
    @views (α |> eltype(A[1, :])) .* A[:, 1] .+ A[:, 2]
end


PDE = SchrodingerPDEPolynomic(boundaries, σ, Nop, Start, T, FieldF)

grid2d = generate_grid(PDE,τ,hx,hy)

A,B,C = generate_ABC(grid2d,6,1.0)

Ker = assembly_metakernel(CPUBackend,get_metadata(grid2d),B,C)

Solver = InitializePDESolver(CPUBackend,get_metadata(grid2d),[Ker;Ker],A)

start = startup_CPU(PDE,grid2d)

full_algorithm(Solver,PDE,grid2d,start)