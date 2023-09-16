#Base Types and functions

include("ComputeBackend/ComputeBackend.jl")
include("AbstractMesh/AbstractMesh.jl")
include("SchrodingerPDE/SchrodingerPDE.jl")
include("SpaceTimeGrid/SpaceTimeGrid.jl")
include("Plot/Plot.jl")
include("Solver/Solver.jl")

#Extended Methods
include("SpaceTimeGrid/methods/index.jl")
include("Plot/methods/index.jl")
include("Solver/methods/index.jl")
