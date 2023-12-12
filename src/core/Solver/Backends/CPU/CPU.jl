export generate_grid, generate_ABC, assembly_metakernel, InitializePDESolver, startup_CPU,
    full_algorithm, update_component!, get_preconditioner_B, solve

function startup_CPU(PDE_Meta::SchrodingerPDE, Grid::SpaceTimeGrid)
    start = zeros(Complex{eltype(get_boundary(PDE_Meta, 1))},
        linear_size(get_metadata(Grid)),
        get_total_components(PDE_Meta))
    for (fun, index) in zip(PDE_Meta.ψ₀, 1:get_total_components(PDE_Meta))
        start[:, index] .= evaluate_indexed_elements(Grid, fun)
    end
    start
end


export plot_inner
function plot_inner(value, mesh::SpaceTimeGrid2D, id::Int = 0)
    indexes = get_indexed_elements(mesh)
    indexes = [x[k] for x in indexes, k in 1:2]

    absv = sum(abs2.(value), dims=2) |> Array
    #absv = abs2.(value[:, 1])

    fig = Figure()
    axe = Axis3(fig[1, 1])
    plot!(axe, indexes[:, 1], indexes[:, 2], absv[:, 1])
    save("pictures/fig$id.png", fig)
    #display(fig)
end

function plot_inner(value, mesh::SpaceTimeGrid1D, id::Int = 0)
    indexes = get_indexed_elements(mesh)
    #indexes = [x[k] for x in indexes, k in 1:2]

    absv = sum(abs2.(value), dims = 2)
    #@show absv
    #absv = abs2.(value[:,1])

    fig = Figure()
    axe = Axis(fig[1, 1])

    plot!(axe, indexes, absv[:, 1])
    save("pictures/fig$id.png", fig)
    #display(fig)
end


@inline function update_component!(Solver::PDESolver2,
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

        broadcast!(*,zcomp,eval_optimized_potential(PDE, current_state, zl, index),zl .+ prev_zl)

        b_temp .-= τ*(opA * zcomp)

        #dst, src
        copy!(zcomp, zl)

        ldiv!(zl, Ker.factorizationB, b_temp)
        

        #Temporary Checks...
        if !all(isfinite,zl)
            throw("Unbounded...")
        end

        zcomp .-= zl

        best_norm = get_square_root_measure(Grid) * norm(zl)
        diff_norm = get_square_root_measure(Grid) * norm(zcomp)



        if diff_norm < a_tol + r_tol * best_norm
            if (l>100)
                plot_inner(current_state,Grid,current_time)
                println("Diff: $(diff_norm) Best: $(diff_norm)")
            end
            update_stats(Stats,0,l)

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


function calculate_N(PDE::SchrodingerPDENonPolynomic, currentMove, updated_component, index::Int)
    nextMove = copy(currentMove)

    nextMove[:,index] .= updated_component

    curr_col = view(currentMove, :, index)
    next_col = updated_component
    diff_ul = abs2.(next_col) .- abs2.(curr_col)

    output = similar(curr_col,eltype(real(curr_col)))

    similar_indexes = abs.(next_col - curr_col) .< 500 * eps(real(eltype(curr_col)))
    non_similar_indexes = .!similar_indexes

    similarities = size(next_col[similar_indexes], 1)
    #@show similarities
    if similarities == 0
        output .= (eval_field(PDE, nextMove) .- eval_field(PDE, currentMove)) ./ (diff_ul)
    elseif similarities == size(curr_col, 1)
        output .= eval_component(PDE, index, currentMove)
    else
        current_similar = view(currentMove, similar_indexes, :)
        current_non_similar = view(currentMove, non_similar_indexes, :)
        diff_non_similar = view(diff_ul, non_similar_indexes, :)
        next_non_similar = view(nextMove, non_similar_indexes, :)

        output[non_similar_indexes] .= (eval_field(PDE, next_non_similar) .-
                                        eval_field(PDE, current_non_similar)) ./
                                       (diff_non_similar)
        output[similar_indexes] .= eval_component(PDE, index, current_similar)
    end

    @show typeof(output)

    output
end

@inline function update_component!(Solver::PDESolver2,
    PDE::SchrodingerPDENonPolynomic,
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

        broadcast!(*,zcomp,calculate_N(PDE, current_state, zl, index),zl .+ prev_zl)

        b_temp .-= τ*(opA * zcomp)

        #dst, src
        copy!(zcomp, zl)

        ldiv!(zl, Ker.factorizationB, b_temp)
        

        #Temporary Checks...
        if !all(isfinite,zl)
            throw("Unbounded...")
        end

        zcomp .-= zl

        best_norm = get_square_root_measure(Grid) * norm(zl)
        diff_norm = get_square_root_measure(Grid) * norm(zcomp)

        #@show best_norm diff_norm r_tol*best_norm+a_tol


        if diff_norm < a_tol + r_tol * best_norm
            if (l>149)
                plot_inner(current_state,Grid,current_time)
                println("Diff: $(diff_norm) Best: $(diff_norm)")
            end
            update_stats(Stats,0,l)

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

function full_algorithmX(start_state, timesteps::Int,
    time_substeps::AbstractArray{T},
    Solver::PDESolver2, Grid::SpaceTimeGrid, PDE::SchrodingerPDE,
    fixed_innerloop_steps::Int,
    r_tol::T,
    a_tol::T, showprogress::Bool, verbose::Int) where {T <: AbstractFloat}

    #Startup Data
    Stats = PDESolverStats(Dictionary{Int64,Int64}(),zeros(Float64,1))
    total_components = get_total_components(PDE)

    σ_f = get_σ(PDE)
    σ_r = reverse(σ_f)

    #Memory Allocations

    current_state = copy(start_state)
    component_place = similar(current_state,
        eltype(current_state),
        size(current_state, 1))

    b_temp = similar(component_place)
    b0_temp = similar(component_place)
    zcomp = similar(component_place)

    #@show size(b0_temp) size(b_temp) size(zcomp) size(current_state) size(component_place) typeof(component_place) typeof(current_state)

    ###############################################
    ##  Main Loop                                ##
    ###############################################
    plot_inner(current_state, Grid, 0)

    start_energy = get_measure(Grid) * sum(abs2.(current_state), dims = 1)
    
    @showprogress "Computing steps" barlen=40 showspeed=true enabled=showprogress for step in 1:timesteps
        @showprogress "Time substep" barlen=40 showspeed=true offset=1 enabled=showprogress for τ in time_substeps
            #Forward 
            for (component_index, σ) in enumerate(σ_f)
                update_component!(Solver,
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
            for (reverse_c_index, σ) in zip(total_components:-1:1,σ_r)
                update_component!(Solver,
                    PDE,
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

        if step % 100 == 0
            @show get_measure(Grid) * sum(abs2.(current_state), dims = 1) - start_energy
            showstats(Stats)
            plot_inner(current_state, Grid, step)
        end
    end
end


function solve(::Type{T},
    ::Type{CPUBackend},
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


    #display(A|>Array)
    #display(D|>Array)


    store_idx = 1
    #Factoring
    @showprogress "Factoring " barlen=40 showspeed=true enabled=showprogress for (σ, βτ) in product(σset,
        TimeMultipliers)
        B = 4im * A + βτ * σ * D
        C = 4im * A - βτ * σ * D

        #display(C |> Array)
        dkeys[store_idx] = (σ, βτ)
        dvalues[store_idx] = MetaKer2(C, lu(B),nothing)
        store_idx += 1
    end

    FactorizationsComponents = Dictionary(dkeys, dvalues)

    #display(FactorizationsComponents)

    PDESolverData = PDESolver2{T, CPUBackend, typeof(Mesh), typeof(A)}(Mesh,
        FactorizationsComponents,
        A)

    #Initialization of solver process...

    state0 = startup_CPU(PDE, Grid)

    TimeCollection = τ * collect(TimeDiscretization)

    @show time_steps
    @show τ
    @show TimeCollection

    full_algorithmX(state0,
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
