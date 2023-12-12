##Approach for access

#Attempt to implement (σ_i,τ_j) -> B,C

@inline function Base.getindex(Solver::PDESolver2, σ::T, τ::T) where {T <: AbstractFloat}
    Solver.Kernel[(σ, τ)]
end

@inline function get_A(Solver::PDESolver2)
    Solver.opA
end

@inline function get_Ker(Solver::PDESolver2, σ::T, τ::T) where {T <: AbstractFloat}
    Solver.Kernel[(σ, τ)]
end

function update_stats(Stats::PDESolverStats, norm, iters)
    if haskey(Stats.IterInfo, iters)
        Stats.IterInfo[iters] += 1
    else
        insert!(Stats.IterInfo, iters, 1)
    end
end

function showstats(Stats::PDESolverStats)
    display(Stats.IterInfo)
end

"Input parsing with some basic checks"
function parse_input_dim_args(::Type{T},
    expected_dims::Int, PDE::SchrodingerPDE;
    kwargs...) where {T <: AbstractFloat}
    h = zeros(T, expected_dims)
    N = zeros(Int64, expected_dims)

    dim_ranges = Array{AbstractRange{T}}(undef, expected_dims)
    time_end = get_end_time(PDE)

    if haskey(kwargs, :τ)
        τ = kwargs[:τ]
        if τ <= 0 || !isa(τ, T)
            throw(DomainError(τ,
                "τ is expected to be a non zero and positive float value with type $(typeof(T))..."))
        end
        time_steps = ceil(Int64, time_end / τ)
    elseif haskey(kwargs, :time_steps)
        time_steps = kwargs[:time_steps]
        if time_steps <= 0 || !isa(time_steps, Int)
            throw(DomainError(time_steps,
                "time_steps is expected to be a non zero and positive integer..."))
        end
        τ = convert(T, time_end / time_steps)

    else
        throw(ArgumentError("You must provide total count of timesteps (time_steps::Int) or time step (τ::$(typeof(T)))..."))
    end

    for (dim, keyh, keyN) in zip([1, 2, 3], [:hx, :hy, :hz], [:Nx, :Ny, :Nz])
        if expected_dims >= dim
            lb, rb = get_boundary(PDE, dim)
            if haskey(kwargs, keyh)
                hval = kwargs[keyh]

                if hval <= 0
                    throw(DomainError(hval,
                        "Solver only works with non zero and positive step $(keyh) in dimension index $(dim)..."))
                end

                drange = lb:hval:rb
                nval = size(drange, 1)

            elseif haskey(kwargs, keyN)
                nval = kwargs[keyN]

                if nval <= 0
                    throw(DomainError(nval,
                        "Solver only works with non zero and positive divisions $(keyN) in dimension with index $(dim)..."))
                end

                drange = range(lb, rb, nval)
                hval = Base.step(drange)

            else
                throw(ArgumentError("Expected $(keyh) or $(keyN). You must define one of them..."))
            end

            N[dim] = nval
            h[dim] = hval
            dim_ranges[dim] = drange
        else
            break
        end
    end

    h, N, dim_ranges, τ, time_steps
end
