##Approach for access

#Attempt to implement (σ_i,τ_j) -> B,C

@inline function Base.getindex(Solver::PDESolver2, σ::T, τ::T) where {T<:AbstractFloat}
    Solver.Kernel[(σ,τ)]
end

@inline function get_A(Solver::PDESolver2)
    Solver.opA
end

@inline function get_Ker(Solver::PDESolver2,σ::T,τ::T) where {T<:AbstractFloat}
    Solver.Kernel[(σ,τ)]
end

function update_stats(Stats::PDESolverStats,norm,iters)
    if haskey(Stats.IterInfo,iters)
        Stats.IterInfo[iters]+=1
    else
        insert!(Stats.IterInfo, iters,1 )
    end
end

function showstats(Stats::PDESolverStats)
    display(Stats.IterInfo)
end