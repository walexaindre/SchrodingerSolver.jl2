using PProf
using SchrodingerSolver
using Profile
function ψ01(x)
    0.2 * (1 - 0.4 * cos(0.5 * x))
end

function ψ02(x, y = 0)
    0.3 * (1 - 0.4 * cos(0.5 * (x + 7.0 * pi / 4.0)))
end

function ψ03(x, y = 0)
    0.5 * (1 - 0.4 * cos(0.5 * x))
end

boundaries = [-4.0 * pi; 4.0 * pi]

σ = [1.0, 1.0, 1.0]

Start = [ψ01, ψ02, ψ03]

T = 100.0

α = 2.0

function N0(prev, current, index::Int)
    p_idx = view(prev, :, index)

    out = 0.5 * (abs2.(p_idx) + abs2.(current))
    if index == 1
        v2 = abs2.(view(prev, :, 2))
        v3 = abs2.(view(prev, :, 3))
        out += (2.0 / 3.0) * v2 + v3
    elseif index == 2
        v1 = abs2.(view(prev, :, 1))
        v3 = abs2.(view(prev, :, 3))
        out += (2.0 / 3.0) * (v1 .+ v3)
    elseif index == 3
        v1 = abs2.(view(prev, :, 1))
        v2 = abs2.(view(prev, :, 2))
        out += v1 .+ (2.0 / 3.0) * v2
    end
    out
end

function FieldF(A)
    0.5 * (A[:, 1] .^ 2 + A[:, 2] .^ 2 + A[:, 3] .^ 2) +
    (2.0 / 3.0) * (A[:, 1] .* A[:, 2] + A[:, 2] .* A[:, 3]) + A[:, 1] .* A[:, 3]
end

PDE = SchrodingerPDEPolynomic(boundaries, σ, N0, Start, T, FieldF)

solve(Float64, CPUBackend, PDE, 2, τ = 0.005, Nx = 2000)
