using SchrodingerSolver
using Test
using Aqua

@testset "SchrodingerSolver.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(SchrodingerSolver)
    end
    include("index.jl")
end
