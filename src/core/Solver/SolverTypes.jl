export SolverMethod,PaulMethod,ADIMethod,PDESolution
abstract type SolverMethod end

abstract type PaulMethod <: SolverMethod end

abstract type ADIMethod <: SolverMethod end

struct PDESolution 
    time_steps::Int #Number of time_steps

end