using PProf
using SchrodingerSolver
using Profile
function ψ01(x, y = 0)
    sqrt(2 |> eltype(x)) * (sech(x + y + 10)) * exp(1im * (x + y) / 4)
end

function ψ02(x, y = 0)
    sqrt(2 |> eltype(x)) * (sech(x - y - 10)) * exp(-1im * (x + y) / 4)
end

boundaries = [-40.0; 40.0]

σ = [1.0,1.0]


Start = [ψ01, ψ02]

T = 40.0

α = 2.0

c = (1+2*α)|>Float64

hx=0.0005
τ = 0.005

function Nop(prevMove,currMove,index::Int)
    #0.5*sum(abs2,currMove,dims=2)+c*view(prevMove,:,3-index)    
    0.5*(abs2.(currMove)+abs2.(view(prevMove,:,index)))+c*abs2.(view(prevMove,:,3-index))   
end

function FieldF(A)
    c₁ = (c |> real(eltype(A)))
    d₁ = (0.5 |> real(eltype(A)))

    B = abs2.(A)
    B₁ = @view B[:, 1]
    B₂ = @view B[:, 2]
    d₁ * (B₁ .^ 2 .+ B₂ .^ 2) .+ c₁ * B₁ .* B₂
end

function ψ1(B)
    A = abs2.(B)
    @views A[:, 1] .+ (c |> eltype(A[1, :])) .* A[:, 2]
end


function ψ2(B)
    A = abs2.(B)
    @views (c |> eltype(A[1, :])) .* A[:, 1] .+ A[:, 2]
end

PDE = SchrodingerPDEPolynomic(boundaries, σ, Nop, Start, T, FieldF)

solve(Float64,CPUBackend,PDE,4,4,5,τ=0.005,Nx=20000)
