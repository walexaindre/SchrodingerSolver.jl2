"Evaluate start conditions..."
function startup_GPU(PDE_Meta::SchrodingerPDE, Grid::SpaceTimeGrid)
    start = zeros(Complex{eltype(PDE_Meta)},
        length(Grid),
        length(PDE_Meta))

    for (fun, index) in zip(PDE_Meta.ψ₀, 1:length(PDE_Meta))
        start[:, index] .= evaluate_indexed_elements(Grid, fun)
    end
    start |> CuArray
end

"Forward - Backward LU z -> Temporary storage, x -> Vector y -> Output placement "
function ldiv_ilu0!(P::CuSparseMatrixCSR, x, y, z)
    ldiv!(z, UnitLowerTriangular(P), x)      # Forward substitution with L
    ldiv!(y, UpperTriangular(P), z)  # Backward substitution with U
    return y
end

"Verify sanity of preconditioner for B with a sample ones vector"
function test_ilu_preconditioner(::Type{T}, Op::LinearOperator) where {T <: AbstractFloat}
    #NaN simple checking of preconditioner (singular ILU)
    if !all(isfinite, Op * CUDA.ones(Complex{T}, size(Op, 1)))
        throw(OverflowError("ILU GPU preconditioner overflows..."))
    end
end

@inline function update_component!(index::Int,
        τ::T,
        σ::T,
        OptionSolver::SchrodingerSolverOptions,
        Memory::BackendMemoryCUDA{T}) where {T <: AbstractFloat}
    PDE = OptionSolver.PDE
    Stats = OptionSolver.stats
    zl = Memory.component_place
    prev_zl = view(Memory.current_state, :, index)
    current_state = Memory.current_state

    b0_temp = Memory.b0_temp
    b_temp = Memory.b_temp
    zcomp = Memory.zcomp
    KrylovSolver = Memory.solver_memory

    a_tol = OptionSolver.a_tol
    r_tol = OptionSolver.r_tol

    Grid = OptionSolver.Solver.Grid

    #dst, src
    copy!(zl, prev_zl)

    opA = get_A(OptionSolver.Solver)
    Ker = get_Ker(OptionSolver.Solver, σ, τ)

    b0_temp .= Ker.C * zl

    success = false

    #@showprogress "Update Component Loop" enabled=showprogress barlen=40 showspeed=true offset=showprogress_offset 
    for l in 1:(OptionSolver.max_iterations)

        #dst, src
        copy!(b_temp, b0_temp)

        broadcast!(*,
            zcomp,
            eval_optimized_potential(PDE, current_state, zl, index),
            zl .+ prev_zl)

        b_temp .+= τ * (opA * zcomp)

        #dst, src
        copy!(zcomp, zl)
        gmres!(KrylovSolver,
            Ker.factorizationB,
            restart = true,
            N = Ker.preconditionerB,
            b_temp,
            atol = a_tol,
            rtol = r_tol)
        copy!(zl, KrylovSolver.x)

        zcomp .-= zl

        #CUDA.@profile norm(zl)

        best_norm = get_square_root_measure(Grid) * norm(zl)
        diff_norm = get_square_root_measure(Grid) * norm(zcomp)
        #@show best_norm diff_norm r_tol*best_norm+a_tol

        if diff_norm < a_tol + r_tol * best_norm
            if (l > 149)
            end
            update_stats(Stats, 0, l)

            success = true
            break
        end
    end

    if !success
        println("Warning... unsucceful convergence after max_iteration: ", max_iterations)

        throw("Error...")
    end

    copy!(prev_zl, zl)

    nothing
end

@inline function update_component!(Solver::PDESolver2,
        KrylovSolver::GmresSolver,
        PDE::SchrodingerPDEPolynomic,
        Grid::SpaceTimeGrid,
        Stats::PDESolverStats,
        b_temp::AbstractArray{Complex{T}},
        b0_temp::AbstractArray{Complex{T}},
        zcomp::AbstractArray{Complex{T}},
        current_state::AbstractArray{Complex{T}},
        component_place::AbstractArray{Complex{T}},
        index::Int,
        current_time::Int,
        τ::T,
        σ::T,
        a_tol::T,
        r_tol::T,
        max_iterations::Int,
        showprogress::Bool,
        showprogress_offset::Int,
        verbose::Int) where {T <: AbstractFloat}

    #assignation

    zl = component_place
    prev_zl = view(current_state, :, index)

    #@show size(zl) size(prev_zl) typeof(zl) typeof(prev_zl)

    #dst, src
    copy!(zl, prev_zl)
    opA = get_A(Solver)
    Ker = get_Ker(Solver, σ, τ)

    b0_temp .= Ker.C * zl

    success = false

    #@showprogress "Update Component Loop" enabled=showprogress barlen=40 showspeed=true offset=showprogress_offset 
    for l in 1:max_iterations

        #dst, src
        copy!(b_temp, b0_temp)

        broadcast!(*,
            zcomp,
            eval_optimized_potential(PDE, current_state, zl, index),
            zl .+ prev_zl)

        b_temp .-= τ * (opA * zcomp)

        #dst, src
        copy!(zcomp, zl)
        gmres!(KrylovSolver,
            Ker.factorizationB,
            restart = true,
            M = Ker.preconditionerB,
            b_temp,
            atol = a_tol,
            rtol = r_tol)
        copy!(zl, KrylovSolver.x)
        #Temporary Checks...
        #if !all(isfinite, zl)
        #    throw("Unbounded...")
        #end

        #display(KrylovSolver.stats)
        zcomp .-= zl

        #CUDA.@profile norm(zl)

        best_norm = get_square_root_measure(Grid) * norm(zl)
        diff_norm = get_square_root_measure(Grid) * norm(zcomp)
        #@show best_norm diff_norm r_tol*best_norm+a_tol

        if diff_norm < a_tol + r_tol * best_norm
            if (l > 149)
                plot_inner(current_state, Grid, current_time)
                println("Diff: $(diff_norm) Best: $(diff_norm)")
            end
            update_stats(Stats, 0, l)

            success = true
            break
        end
    end

    if !success
        println("Warning... unsucceful convergence after max_iteration: ", max_iterations)

        throw("Error...")
    end

    copy!(prev_zl, zl)

    nothing
end

function full_algorithmY(start_state, timesteps::Int,
        time_substeps::AbstractArray{T},
        Solver::PDESolver2, Grid::SpaceTimeGrid, PDE::SchrodingerPDE,
        fixed_innerloop_steps::Int,
        r_tol::T,
        a_tol::T, showprogress::Bool, verbose::Int) where {T <: AbstractFloat}

    #Startup Data
    Stats = PDESolverStats(Dictionary{Int64, Int64}(), zeros(Float64, 1))
    total_components = length(PDE)

    σ_f = get_σ(PDE)
    σ_r = reverse(σ_f)

    n = length(Grid)
    #Memory Allocations

    current_state = copy(start_state)
    component_place = similar(current_state,
        eltype(current_state),
        n)

    b_temp = similar(component_place)
    b0_temp = similar(component_place)
    zcomp = similar(component_place)

    KrylovSolver = GmresSolver(n, n, 20, CuVector{Complex{T}})

    #@show size(b0_temp) size(b_temp) size(zcomp) size(current_state) size(component_place) typeof(component_place) typeof(current_state)

    ###############################################
    ##  Main Loop                                ##
    ###############################################
    plot_inner(current_state, Grid, 0)

    start_energy = get_measure(Grid) * sum(abs2.(current_state), dims = 1)

    @showprogress "GPU Computing steps" barlen=40 showspeed=true enabled=showprogress for step in 1:timesteps
        @showprogress "Time substep" barlen=40 showspeed=true offset=1 enabled=showprogress for τ in time_substeps
            #Forward 
            #io = IOBuffer()
            #CUDA.@profile io=io 
            begin
                for (component_index, σ) in enumerate(σ_f)
                    update_component!(Solver,
                        KrylovSolver,
                        PDE,
                        Grid,
                        Stats,
                        b_temp,
                        b0_temp,
                        zcomp,
                        current_state,
                        component_place,
                        component_index,
                        step,
                        τ,
                        σ,
                        a_tol,
                        r_tol,
                        150,
                        showprogress,
                        2,
                        verbose)
                end

                #Backward
                for (reverse_c_index, σ) in zip(total_components:-1:1, σ_r)
                    update_component!(Solver,
                        KrylovSolver, PDE,
                        Grid,
                        Stats,
                        b_temp,
                        b0_temp,
                        zcomp,
                        current_state,
                        component_place,
                        reverse_c_index,
                        step,
                        τ,
                        σ,
                        a_tol,
                        r_tol,
                        150,
                        showprogress,
                        3,
                        verbose)
                end
            end

            #filebuff = open("profile.txt", "w+")
            #write(filebuff, take!(io))

            #close(filebuff)
            #close(io)
            #throw("io end")
        end

        if step % 100 == 0
            @show get_measure(Grid) * sum(abs2.(current_state), dims = 1) - start_energy
            showstats(Stats)
            plot_inner(current_state, Grid, step)
        end
    end
end

function initialize(::Type{T},
        ::Type{GPUBackend},
        PDE::SchrodingerPDE,
        space_order::Int = 2,
        time_order::Int = 2,
        time_composition_substeps::Int = 1,
        time_composition_index::Int = 1;
        fixed_innerloop_steps::Int = 0,
        r_tol::T = 700 * eps(T),
        a_tol::T = 700 * eps(T),
        showprogress::Bool = true,
        verbose::Int = 0,
        kwargs...) where {T <: AbstractFloat}
    Grid, TimeCollection, TimeMultipliers, time_steps, σset = parse_input_parameters(T,
        PDE,
        time_order,
        time_composition_substeps,
        time_composition_index; kwargs...)

    Mesh = get_metadata(Grid)
    A = get_sparsematrix_A(T, Mesh, space_order)
    D = get_sparsematrix_D(Grid, space_order)

    #We need to generate B, C for all time multipliers....
    dkeys = Array{Tuple{T, T}}(undef, length(σset) * length(TimeMultipliers))
    dvalues = Array{MetaKer2}(undef, length(σset) * length(TimeMultipliers))

    store_idx = 1
    #Preconditioner and Calculations

    @showprogress "Preconditioner start and checks" barlen=40 showspeed=true enabled=showprogress for (σ, βτ) in product(σset,
        TimeMultipliers)
        B = (4im * A + βτ * σ * D) |> CuSparseMatrixCSR

        opM = sparse(drop(B, Mesh)..., fmt = :csr)

        C = (4im * A - βτ * σ * D) |> CuSparseMatrixCSR

        dkeys[store_idx] = (σ, βτ)
        dvalues[store_idx] = MetaKer2(C, B, opM)
        store_idx += 1
    end

    A = convert(SparseMatrixCSC{ComplexF64}, A) |> CuSparseMatrixCSR

    D = convert(SparseMatrixCSC{ComplexF64}, D) |> CuSparseMatrixCSR
  
    FactorizationsComponents = Dictionary(dkeys, dvalues)

    opA = sparse(drop(A, Mesh)..., fmt = :csr)
    PDESolverData = PDESolver3{T, GPUBackend, typeof(Grid), typeof(A)}(Grid,
        FactorizationsComponents,
        D,
        A,
        opA)


    Memory = initialize_memory(T, length(PDE), length(Mesh), 20)


    Memory.current_state .= startup_GPU(PDE, Grid)

    "is_initialized = true,
    r_tol = r_tol,
    a_tol = a_tol,
    fixed_innerloop_steps = fixed_innerloop_steps,
    max_iterations = max_iterations,
    showprogress = showprogress,
    verbose = verbose,
    stats=initialize_statics(),
    PDE=PDE,
    Solver = PDESolverData,
    start_power = 0
    start_energy = [0,0]
    time_collection=TimeCollection,
    compute_backend = GPUBackend,
    data_type = T"

    start_power = system_power(Grid, Memory)
    start_energy = T(0)

    @show start_power

    Opt = SchrodingerSolverOptions{T, GPUBackend}(true,
        r_tol,
        a_tol,
        fixed_innerloop_steps,
        200,
        showprogress,
        verbose,
        initialize_statics(),
        PDE,
        PDESolverData,
        start_power,
        start_energy,
        TimeCollection,
        GPUBackend,
        T)

    Opt.start_energy = system_energy(Opt, Memory)

    Opt, Memory
end

"GPU Step Solver"
function step!(OptionSolver::SchrodingerSolverOptions, Memory::BackendMemoryCUDA)
    σ_f = get_σ(OptionSolver.PDE)
    σ_r = reverse(σ_f)

    for τ in OptionSolver.time_collection
        begin
            #Forward
            for (component_index, σ) in enumerate(σ_f)
                update_component!(component_index, τ, σ, OptionSolver, Memory)
            end

            #Backward
            for (reverse_c_index, σ) in zip(length(σ_r):-1:1, σ_r)
                update_component!(reverse_c_index, τ, σ, OptionSolver, Memory)
            end
        end
    end
end

"GPU Solver..."
function solve(::Type{T},
        ::Type{GPUBackend},
        PDE::SchrodingerPDE,
        space_order::Int = 2,
        time_order::Int = 2,
        time_composition_substeps::Int = 1,
        time_composition_index::Int = 1;
        fixed_innerloop_steps::Int = 0,
        r_tol::T = 700 * eps(T),
        a_tol::T = 700 * eps(T),
        showprogress::Bool = true,
        verbose::Int = 0,
        kwargs...) where {T <: AbstractFloat}

    #Dimensionality
    expected_dims = ndims(PDE)

    #Validations and startup calculations...
    h, N, dim_ranges, τ, time_steps = parse_input_dim_args(T,
        expected_dims,
        PDE;
        kwargs...)

    @show h N

    #Abstract Mesh with data...
    Grid = CreateGrid(dim_ranges..., τ, h..., N...)
    Mesh = get_metadata(Grid)

    #Initialization of Time Discretization

    TimeDiscretization = get_time_discretization(T,
        time_order,
        time_composition_substeps,
        time_composition_index)

    #Initialization of Space Discretization

    TimeMultipliers = τ * get_coefficients(TimeDiscretization)

    σset = Set(get_σ(PDE))

    #We need to generate B, C for all time multipliers....
    dkeys = Array{Tuple{T, T}}(undef, length(σset) * length(TimeMultipliers))
    dvalues = Array{MetaKer2}(undef, length(σset) * length(TimeMultipliers))

    A = get_sparsematrix_A(T, Mesh, space_order)
    D = get_sparsematrix_D(Grid, space_order)
    nsz = size(A, 1)

    store_idx = 1
    #Preconditioner and Calculations

    @showprogress "Preconditioner start and checks" barlen=40 showspeed=true enabled=showprogress for (σ, βτ) in product(σset,
        TimeMultipliers)
        B = (4im * A + βτ * σ * D) |> CuSparseMatrixCSR

        opM = sparse(drop(B, Mesh)..., fmt = :csr)
        @show "B sz" size(B) "opM sz" size(opM)
        C = (4im * A - βτ * σ * D) |> CuSparseMatrixCSR

        dkeys[store_idx] = (σ, βτ)
        dvalues[store_idx] = MetaKer2(C, B, opM)
        store_idx += 1
    end
    A = convert(SparseMatrixCSC{ComplexF64}, A) |> CuSparseMatrixCSR
    FactorizationsComponents = Dictionary(dkeys, dvalues)

    #display(FactorizationsComponents)

    PDESolverData = PDESolver2{T, GPUBackend, typeof(Mesh), typeof(A)}(Mesh,
        FactorizationsComponents,
        A)

    #Initialization of solver process...

    state0 = startup_GPU(PDE, Grid)

    TimeCollection = τ * collect(TimeDiscretization)

    @show time_steps
    @show τ
    @show TimeCollection

    full_algorithmY(state0,
        time_steps,
        TimeCollection,
        PDESolverData,
        Grid,
        PDE,
        fixed_innerloop_steps,
        r_tol,
        a_tol,
        showprogress,
        verbose)
end