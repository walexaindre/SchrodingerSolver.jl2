using SchrodingerSolver
using Test

@testset "SchrodingerSolver.jl" begin
    # Write your tests here.
    @test  compute_backend=="CUDA"
end
