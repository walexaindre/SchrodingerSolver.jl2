@inline function system_power(Grid::SpaceTimeGrid, Memory::BackendMemory)
    get_measure(Grid) * vec(sum(abs2.(Memory.current_state), dims = 1))
end

@inline function system_energy(Opt::SchrodingerSolverOptions{T},
        Memory::BackendMemory) where {T <: AbstractFloat}
    preA = get_preconditionerA(Opt.Solver)
    D = get_D(Opt.Solver)
    A = get_A(Opt.Solver)

    components = Memory.current_state

    energy = T(0)

    for (idx, σ) in enumerate(get_σ(Opt.PDE))
        comp = components[:, idx]
        b = D * comp

        gmres!(Memory.solver_memory,
            A,
            b,
            restart = true,
            N = preA,
            atol = Opt.a_tol,
            rtol = Opt.r_tol)

        energy -= σ * dot(comp, Memory.solver_memory.x)
    end

    v2 = sum(eval_field(Opt.PDE, components), dims = 1) |> Array
    v2 = v2[1]
    energy += v2

    T(0.5) * get_measure(Opt.Solver.Grid) * real(energy)
end

function Error(Opt::SchrodingerSolverOptions, Memory::BackendMemory)
    abs.((system_power(Opt.Solver.Grid, Memory)|>typeof(Opt.start_power)) - Opt.start_power),
    abs(system_energy(Opt, Memory) - Opt.start_energy)
end