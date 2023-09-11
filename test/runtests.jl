using ShrodingerSolver
using Test

@testset "ShrodingerSolver.jl" begin
    # Write your tests here.
    @test  compute_backend=="CUDA"
end
