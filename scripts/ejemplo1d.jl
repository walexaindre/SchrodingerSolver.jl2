using SchrodingerSolver


#=

=#

function ψ01(x, y = 0)
    sqrt(2 |> eltype(x)) * (sech(x + y + 10)) * exp(1im * (x + y) / 4)
end

function ψ02(x, y = 0)
    sqrt(2 |> eltype(x)) * (sech(x - y - 10)) * exp(-1im * (x + y) / 4)
end

boundaries = [-40.0; 40.0]

σ = [1.0; 1.0]


Start = [ψ01, ψ02]

T = 40.0

α = 1.0

c = (1+2*α)|>Float64

hx=0.25
hy=0.5
τ = 0.005

function Nop(prevMove,currMove,index::Int)
    #0.5*sum(abs2,currMove,dims=2)+c*view(prevMove,:,3-index)    
    0.5*(abs2.(view(currMove,:,index))+abs2.(view(prevMove,:,index)))+c*abs2.(view(prevMove,:,3-index))   
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

grid1d = generate_grid(PDE,τ,hx)

cpu = startup_CPU(PDE,grid1d)

kerb = generate_kernel_B(Float64,get_metadata(grid1d),τ,T,hx)

kerc = generate_kernel_C(Float64,get_metadata(grid1d),τ,T,hx)

preb = generate_preconditioner_B(kerb)

metaker = assembly_metakernel(CPUBackend,get_metadata(grid1d),kerb,kerc,preb)



solver = PDESolver{CPUBackend,MetaMesh1D}(get_metadata(grid1d),[metaker,metaker])

full_algorithm(solver,PDE,grid1d,cpu)