using PProf
using SchrodingerSolver
using Profile
function ψ01(x)
    sqrt(2.0) * sech(x + 20) * exp(x * 1im / 4.0)
end

function ψ02(x)
    sqrt(2.0) * sech(x + 10) * exp(x * 1im / 4.0)
end

function ψ03(x)
    sqrt(2.0) * sech(x - 10) * exp(-x * 1im / 4.0)
end

function ψ04(x)
    sqrt(2.0) * sech(x - 20) * exp(-x * 1im / 4.0)
end

boundaries = [-40.0; 40.0]

σ = [1.0, 1.0, 1.0, 1.0]

Start = [ψ01, ψ02, ψ03, ψ04]

T = 50.0

α = 1.0
κ = 1.0
e_c = 3.0

function N0(prev, current, index::Int)
    p_idx = view(prev, :, index)
    k_idx = view(prev, :, 5 - index)
    out = (α/2) * (abs2.(p_idx) .+ abs2.(current)) .+ κ * k_idx

    if index == 1 || index == 4
        v2 = abs2.(view(prev, :, 2))
        v3 = abs2.(view(prev, :, 3))
        out .+= e_c * (v2 .+ v3)
    elseif index == 2 || index == 3
        v1 = abs2.(view(prev, :, 1))
        v4 = abs2.(view(prev, :, 4))
        out .+= e_c * (v1 .+ v4)
    else
        throw("Bound error...")
    end
    out
end

function FieldF(A)
    x = abs2.(view(A, :, 1))
    y = abs2.(view(A, :, 2))
    z = abs2.(view(A, :, 3))
    w = abs2.(view(A, :, 4))
    @. (α / 2) * (x^2 + y^2 + z^2 + w^2) + e_c * (x * y + x * z + y * w + z * w) +
       κ * (x * w + y * z)
end

PDE = SchrodingerPDEPolynomic(boundaries, σ, N0, Start, T, FieldF)

solve(Float64, CPUBackend, PDE, 2, τ = 0.00625, hx = 0.1)
