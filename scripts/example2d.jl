using SchrodingerSolver


#=

=#

function ψ01(x, y = 0)
    sqrt(2 |> eltype(x)) * (sech(x + y + 10)) * exp(1im * (x + y) / 4)
end

function ψ02(x, y = 0)
    sqrt(2 |> eltype(x)) * (sech(x - y - 10)) * exp(-1im * (x + y) / 4)
end

boundaries = [-40.0 -40.0; 40.0 40.0]

σ = [1.0; 1.0]


Start = [ψ01, ψ02]

T = 40.0

α = 1.0

c = (1+2*α)

hx=0.5
hy=0.5
τ = 0.05
function Nop(prevMove,currMove,index::Int)
    0.5*sum(abs2,currMove,dims=2)+c*view(prevMove,:,3-index)    
end

function FieldF(A)
    c₁ = (c |> real(eltype(A)))
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