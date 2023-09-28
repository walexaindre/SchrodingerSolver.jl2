export generate_grid, generate_ABC, assembly_metakernel, InitializePDESolver, startup_CPU,
    full_algorithm, update_component!, get_preconditioner_B, solve

function generate_grid(PDE::SchrodingerPDE, τ::T, hx::T) where {T <: AbstractFloat}
    function get_space(min_space, max_space, dist)
        collect(min_space:dist:max_space)
    end

    dim = ndims(PDE)

    if dim == 1 && hx != 0
        l, r = get_boundary(PDE, 1)
        space = get_space(l, r - hx, hx)
        mesh = MetaMesh1D(size(space, 1))
    else
        throw(DimensionMismatch("Dimension or invalid step size dim: $(dim) expected dim: 1 step: $(hx)"))
    end

    SpaceTimeGrid1D(space, mesh, τ, hx)
end

function generate_grid(PDE::SchrodingerPDE, τ::T, hx::T, hy::T) where {T <: AbstractFloat}
    function get_space(min_space, max_space, dist)
        collect(min_space:dist:max_space)
    end

    dim = ndims(PDE)

    if dim == 2 && hx != 0 && hy != 0
        lx, rx = get_boundary(PDE, 1)
        ly, ry = get_boundary(PDE, 2)

        spacex = lx:hx:rx
        spacey = ly:hy:ry

        mesh = MetaMesh2D(size(spacex, 1), size(spacey, 1))
    else
        throw(DimensionMismatch("Dimension or invalid step size dim: $(dim) expected dim: 2 steps: x-> $(hx) y-> $(hy)"))
    end

    SpaceTimeGrid2D(spacex, spacey, mesh, τ, hx, hy)
end

function generate_grid(PDE::SchrodingerPDE,
    τ::T,
    hx::T,
    hy::T,
    hz::T) where {T <: AbstractFloat}
    function get_space(min_space, max_space, dist)
        collect(min_space:dist:max_space)
    end

    dim = ndims(PDE)

    if dim == 3 && hx != 0 && hy != 0 && hz != 0
        lx, rx = get_boundary(PDE, 1)
        ly, ry = get_boundary(PDE, 2)
        lz, rz = get_boundary(PDE, 3)

        spacex = get_space(lx, rx - hx, hx)
        spacey = get_space(ly, ry - hy, hy)
        spacez = get_space(lz, rz - hz, hz)

        mesh = MetaMesh3D(size(spacex, 1), size(spacey, 1), size(spacez, 1))
    else
        throw(DimensionMismatch("Dimension or invalid step size dim: $(dim) expected dim: 3, steps: x-> $(hx) y-> $(hy)"))
    end

    SpaceTimeGrid3D(spacex, spacey, spacez, mesh, τ, hx, hy, hz)
end

function generate_ABC(Grid::SpaceTimeGrid, order::Int, σ::T) where {T <: AbstractFloat}
    A = get_sparsematrix_A(eltype(Grid.hx), get_metadata(Grid), order)
    D = get_sparsematrix_D(Grid, order)

    D .*= Grid.τ * σ
    D = convert(SparseMatrixCSC{Complex{eltype(D)}}, D)
    A = convert(SparseMatrixCSC{Complex{eltype(A)}}, A)
    Aout = deepcopy(A)
    A .*= 4im
    Aout, A + D, A - D
end

function assembly_metakernel(::Type{CPUBackend},
    Mesh::MetaMesh,
    B::SparseMatrixCSC,
    C::SparseMatrixCSC)
    lp = LinearProblem(B, ones(eltype(B), size(B, 1)))

    linsolve = init(lp, UMFPACKFactorization())
    #cond(B,1)
    #@show linsolve
    MetaKernel(B, C, linsolve)
end

function InitializePDESolver(::Type{Backend},
    Grid::MetaMesh,
    Kernels,
    OperatorA::OperatorOrMatrix = I) where {Backend <: AbstractBackend}
    PDESolver{Backend, typeof(Grid), typeof(OperatorA)}(Grid, Kernels, OperatorA)
end

function startup_CPU(PDE_Meta::SchrodingerPDE, Grid::SpaceTimeGrid)
    start = zeros(Complex{eltype(get_boundary(PDE_Meta, 1))},
        linear_size(get_metadata(Grid)),
        get_total_components(PDE_Meta))
    for (fun, index) in zip(PDE_Meta.ψ₀, 1:get_total_components(PDE_Meta))
        start[:, index] .= evaluate_indexed_elements(Grid, fun)
    end
    start
end

function calculate_N(PDE::SchrodingerPDENonPolynomic, currentMove, nextMove, index::Int)
    curr_col = view(currentMove, :, index)
    next_col = view(nextMove, :, index)
    diff_ul = abs2.(next_col) .- abs2.(curr_col)

    output = similar(curr_col)

    similar_indexes = abs.(next_col - curr_col) .< 500 * eps(real(eltype(curr_col)))
    non_similar_indexes = .!similar_indexes

    similarities = size(next_col[similar_indexes], 1)
    #@show similarities
    if similarities == 0
        output .= (eval_field(PDE, nextMove) .- eval_field(PDE, currentMove)) ./ (diff_ul)
    elseif similarities == size(similarities, 1)
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
    output
end

function get_operator_A(Solver::PDESolver)
    return Solver.opA
end

function get_kernel_B(Solver::PDESolver, index::Int)
    Solver.Kernel[index].linearOperatorB
end

function get_kernel_C(Solver::PDESolver, index::Int)
    Solver.Kernel[index].linearOperatorC
end

function get_preconditioner_B(Solver::PDESolver, index::Int)
    Solver.Kernel[index].preconditionerB
end

function update_component!(PDEsln::PDESolver,
    PDE::SchrodingerPDENonPolynomic,
    Grid::SpaceTimeGrid,
    currentMove,
    nextMove,
    index::Int,
    max_iterations::Int = 500,
    error::AbstractFloat = 7e-14)
    #opB = get_kernel_B(PDEsln, index)
    opA = get_operator_A(PDEsln)
    opC = get_kernel_C(PDEsln, index)
    preB = get_preconditioner_B(PDEsln, index)

    zl = view(nextMove, :, index)
    prev_zl = view(currentMove, :, index)

    b0 = opC * prev_zl
    zcomp = similar(b0)
    #b = similar(b0)
    τ = get_τ(Grid)

    success = false

    for l in 1:max_iterations
        zcomp .= zl

        preB.b = b0 .-
                 τ .* opA *
                 (calculate_N(PDE, currentMove, nextMove, index) .* (zl .+ prev_zl))

        sol = LinearSolve.solve!(preB)
        zl .= sol.u

        #gmres!(solver,opB,b)

        #zl.=solver.x
        zcomp .-= zl

        if (l > 1500)
            @show get_square_root_measure(Grid)
            @show norm(zcomp)
            @show get_square_root_measure(Grid) * norm(zcomp)
        end

        if get_square_root_measure(Grid) * norm(zcomp) < error
            success = true
            break
        end
    end

    #@show norm(currentMove[:,1]-nextMove[:,1]) norm(currentMove[:,2]-nextMove[:,2])
    if !success
        println("Warning... unsucceful convergence after max_iteration: ", max_iterations)
        throw("Error...")
    end

    prev_zl .= zl
    nothing
end

function update_component!(PDEsln::PDESolver,
    PDE::SchrodingerPDEPolynomic,
    Grid::SpaceTimeGrid,
    currentMove,
    nextMove,
    index::Int,
    max_iterations::Int = 500,
    error::AbstractFloat = 7e-14)
    opA = get_operator_A(PDEsln)
    #opB = get_kernel_B(PDEsln, index)
    opC = get_kernel_C(PDEsln, index)
    preB = get_preconditioner_B(PDEsln, index)

    zl = view(nextMove, :, index)
    prev_zl = view(currentMove, :, index)

    b0 = opC * prev_zl
    zcomp = similar(b0)
    #b = similar(b0)
    τ = get_τ(Grid)

    success = false

    for l in 1:max_iterations
        zcomp .= zl

        preB.b = b0 .-
                 τ .* opA *
                 (eval_optimized_potential(PDE, currentMove, nextMove, index) .*
                  (zl .+ prev_zl))

        sol = LinearSolve.solve!(preB)
        zl .= sol.u

        #gmres!(solver,opB,b)

        #zl.=solver.x
        zcomp .-= zl

        if (l > 1500)
            @show get_square_root_measure(Grid)
            @show norm(zcomp)
            @show get_square_root_measure(Grid) * norm(zcomp)
        end

        if get_square_root_measure(Grid) * norm(zcomp) < error
            if l > 50
                @show l
            end
            success = true
            break
        end
    end

    #@show norm(currentMove[:,1]-nextMove[:,1]) norm(currentMove[:,2]-nextMove[:,2])
    if !success
        println("Warning... unsucceful convergence after max_iteration: ", max_iterations)
        throw("Error...")
    end

    prev_zl .= zl
    nothing
end

export plot_inner
function plot_inner(value, mesh::SpaceTimeGrid2D, id::Int = 0)
    indexes = get_indexed_elements(mesh)
    indexes = [x[k] for x in indexes, k in 1:2]

    #absv = sum(abs2.(value), dims=2) |> Array
    absv = abs2.(value[:, 1])

    fig = Figure()
    axe = Axis3(fig[1, 1])
    plot!(axe, indexes[:, 1], indexes[:, 2], absv[:, 1])
    save("pictures/fig$id.png", fig)
    display(fig)
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
    display(fig)
end

function full_algorithm(Solver::PDESolver, PDE::SchrodingerPDE, Grid::SpaceTimeGrid, Ψⁿ)
    time_steps = Int64(get_time_steps(PDE, Grid))

    total_components = get_total_components(PDE)

    #Main idea, no more allocations more than this!

    currentMove = copy(Ψⁿ)
    nextMove = copy(Ψⁿ)

    #element_count = linear_size(get_metadata(grid))
    #solver = GmresSolver(element_count,element_count,20,typeof(currentMove[:,1]))

    #End of allocations..
    @show size(currentMove)
    @show time_steps
    plot_inner(currentMove, Grid, 0)
    start_energy = get_measure(Grid) * sum(abs2.(currentMove), dims = 1)
    @showprogress "Computing" showspeed=true for time_step in 1:time_steps

        #Forward
        for component_index in 1:total_components
            update_component!(Solver, PDE, Grid, currentMove, nextMove, component_index)
        end

        #@show "pass2"
        #Backward
        for reverse_c_index in total_components:-1:1
            update_component!(Solver, PDE, Grid, currentMove, nextMove, reverse_c_index)
        end

        #display(currentMove)
        if time_step % 100 == 0
            plot_inner(currentMove, Grid, time_step)
            @show get_measure(Grid) * sum(abs2.(currentMove), dims = 1) - start_energy
            @show time_step
        end
    end

    return nextMove
end

function check_consistency()
end

function calculate_dim_h(a::T, b::T, N::Int) where {T <: AbstractFloat}
    (b - a) / (N - 1)
end

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

                drange = lb:hval:rb
                nval = size(drange, 1)

                if hval <= 0
                    throw(DomainError(hval,
                        "Solver only works with non zero and positive step $(keyh) in dimension index $(dim)..."))
                end
            elseif haskey(kwargs, keyN)
                nval = kwargs[keyN]
                hval = calculate_dim_h(lb, rb, nval)

                if nval <= 0
                    throw(DomainError(nval,
                        "Solver only works with non zero and positive divisions $(keyN) in dimension with index $(dim)..."))
                end
            else
                throw(ArgumentError("Expected $(keyh) or $(keyN). You must define one of them..."))
            end

            N[dim] = nval
            h[dim] = hval
            dim_ranges[dim] = lb:hval:rb
        else
            break
        end
    end

    h, N, dim_ranges, τ, time_steps
end

function solve(::Type{T},
    ::Type{CPUBackend},
    PDE::SchrodingerPDE,
    space_order::Int = 2,
    time_order::Int = 2;
    time_composition_substeps::Int = 1,
    time_composition_index::Int = 1,
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

    store_idx = 1
    for (σ, βτ) in product(σset, TimeMultipliers)
        B = 4im * A - βτ * σ * D
        C = 4im * A + βτ * σ * D
        dkeys[store_idx] = (σ, βτ)
        dvalues[store_idx] = MetaKer2(C, lu(B))
    end

    FactorizationsComponents = Dictionary(dkeys, dvalues)

    PDESolverData = PDESolver2{T, CPUBackend, typeof(Mesh), typeof(A)}(Mesh,
        FactorizationsComponents,
        A)

    #Initialization of solver process...

    state0 = startup_CPU(PDE, Grid)

    sz = 1
end

function full_algorithmX(start_state,
    Solver::PDESolver2,
    fixed_innerloop_steps::Int,
    r_tol::T,
    a_tol::T, showprogress::Bool, verbose::Int) where {T <: AbstractFloat}

    #Memory Allocations

end

#=

function full_algorithm(Solver::PDESolver, PDE::SchrodingerPDE, Grid::SpaceTimeGrid, Ψⁿ)
    time_steps = Int64(get_time_steps(PDE, Grid))

    total_components = get_total_components(PDE)

    #Main idea, no more allocations more than this!

    currentMove = copy(Ψⁿ)
    nextMove = copy(Ψⁿ)

    #element_count = linear_size(get_metadata(grid))
    #solver = GmresSolver(element_count,element_count,20,typeof(currentMove[:,1]))

    #End of allocations..
    @show size(currentMove)
    @show time_steps
    plot_inner(currentMove, Grid, 0)
    start_energy = get_measure(Grid) * sum(abs2.(currentMove), dims = 1)
    @showprogress "Computing" showspeed=true for time_step in 1:time_steps

        #Forward
        for component_index in 1:total_components
            update_component!(Solver, PDE, Grid, currentMove, nextMove, component_index)
        end

        #@show "pass2"
        #Backward
        for reverse_c_index in total_components:-1:1
            update_component!(Solver, PDE, Grid, currentMove, nextMove, reverse_c_index)
        end

        #display(currentMove)
        if time_step % 100 == 0
            plot_inner(currentMove, Grid, time_step)
            @show get_measure(Grid) * sum(abs2.(currentMove), dims = 1) - start_energy
            @show time_step
        end
    end

    return nextMove
end
=#