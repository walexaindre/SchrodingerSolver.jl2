using SchrodingerSolver


#=

=#

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

ψ02 = ψ01

boundaries = [-20.0 -20.0; 20.0 20.0]

σ = [1.0; 1.0]


Start = [ψ01, ψ02]

T = 0.1

α = 1.0

hx=0.05
hy=0.05

τ = 0.00000625
function Nop(prevMove,currMove,index::Int)
    0.5*sum(abs2,currMove,dims=2)+α*view(prevMove,:,3-index) 
end

function FieldF(A)
    c₁ = (α |> real(eltype(A)))
    d₁ = (0.5 |> real(eltype(A)))

    B = abs2.(A)
    B₁ = @view B[:, 1]
    B₂ = @view B[:, 2]
    d₁ * (B₁ .^ 2 .+ B₂ .^ 2) .+ c₁ * B₁ .* B₂
end



PDE = SchrodingerPDEPolynomic(boundaries, σ, Nop, Start, T, FieldF)

grid2d = generate_grid(PDE,τ,hx,hy)

cpu = startup_CPU(PDE,grid2d)

kerb = generate_kernel_B(Float64,get_metadata(grid2d),τ,T,hx,hy)

kerc = generate_kernel_C(Float64,get_metadata(grid2d),τ,T,hx,hy)

preb = generate_preconditioner_B(kerb)

metaker = assembly_metakernel(CPUBackend,get_metadata(grid2d),kerb,kerc,preb)



solver = PDESolver{CPUBackend,MetaMesh2D}(get_metadata(grid2d),[metaker,metaker])

full_algorithm(solver,PDE,grid2d,cpu)