include("types.jl")
export generate_kernel_B,
    generate_kernel_C,
    generate_preconditioner_B,
    generate_grid,
    startup_CPU,
    update_component!,
    full_algorithm


function get_kernel_B(Solver::PDESolver, index::Int)
    Solver.Kernel[index].linearOperatorB
end

function get_kernel_C(Solver::PDESolver, index::Int)
    Solver.Kernel[index].linearOperatorC
end

function get_preconditioner_B(Solver::PDESolver, index::Int)
    Solver.Kernel[index].preconditionerB
end



function generate_grid(PDE::SchrodingerPDEBase, τ::AbstractFloat, hx::AbstractFloat)
    τ, hx = promote(τ, hx)
    Type = eltype(get_boundary(PDE, 1))
    if !isa(τ, Type)
        throw("Type Error, your data must have a unique type $(Type)")
    end

    function get_space(min_space, max_space, dist)
        collect(min_space:dist:max_space)
    end

    dim = get_space_dimension(PDE)

    if dim == 1 && hx !== 0
        l, r = get_boundary(PDE, 1)
        space = get_space(l, r, hx)
        mesh = MetaMesh1D(size(space, 1))
        return SpaceTimeGrid1D(space, mesh, τ, hx)
    end

    throw("Dimension or step error...")

end

function generate_grid(
    PDE::SchrodingerPDEBase,
    τ::AbstractFloat,
    hx::AbstractFloat,
    hy::AbstractFloat,
)
    τ, hx, hy = promote(τ, hx, hy)
    Type = eltype(get_boundary(PDE, 1))
    if !isa(τ, Type)
        throw("Type Error, your data must have a unique type $(Type)")
    end

    function get_space(min_space, max_space, dist)
        collect(min_space:dist:max_space)
    end

    dim = get_space_dimension(PDE)

    if dim == 2 && hx !== 0 && hy !== 0

        l1, r1 = get_boundary(PDE, 1)
        l2, r2 = get_boundary(PDE, 2)

        spacex = get_space(l1, r1, hx)
        spacey = get_space(l2, r2, hy)

        mesh = MetaMesh2D(size(spacex, 1), size(spacey, 1))

        return SpaceTimeGrid2D(spacex, spacey, mesh, τ, hx, hy)
    end

    throw("Dimension or step error...")

end

function generate_grid(
    PDE::SchrodingerPDEBase,
    τ::AbstractFloat,
    hx::AbstractFloat,
    hy::AbstractFloat,
    hz::AbstractFloat,
)
    τ, hx, hy, hz = promote(τ, hx, hy, hz)
    Type = eltype(get_boundary(PDE, 1))
    if !isa(τ, Type)
        throw("Type Error, your data must have a unique type $(Type)")
    end

    function get_space(min_space, max_space, dist)
        collect(min_space:dist:max_space)
    end

    dim = get_space_dimension(PDE)

    if dim == 3 && hx !== 0 && hy !== 0 && hz !== 0
        l1, r1 = get_boundary(PDE, 1)
        l2, r2 = get_boundary(PDE, 2)
        l3, r3 = get_boundary(PDE, 3)

        spacex = get_space(l1, r1, hx)
        spacey = get_space(l2, r2, hy)
        spacez = get_space(l3, r3, hz)

        mesh = MetaMesh3D(size(spacex, 1), size(spacey, 1), size(spacez, 1))

        return SpaceTimeGrid3D(spacex, spacey, spacez, mesh, τ, hx, hy, hz)
    end

    throw("Dimension or step error...")
end



function update_component!(
    PDEsln::PDESolver,
    PDE::SchrodingerPDEPolynomic,
    Grid::SpaceTimeGrid,
    currentMove,
    nextMove,
    index::Int,
    max_iterations::Int = 500,
    error::AbstractFloat = 1e-8,
)
    #opB = get_kernel_B(PDEsln, index)
    opC = get_kernel_C(PDEsln, index)
    preB = get_preconditioner_B(PDEsln, index)

    zl = view(nextMove, :, index)
    prev_zl = view(currentMove, :, index)

    b0 = opC * prev_zl
    zcomp = similar(b0)
    b = similar(b0)
    τ = get_τ(Grid)

    success = false

    for l = 1:max_iterations

        zcomp .= zl

        b .=
            b0 .-
            τ .*
            (eval_optimized_potential(PDE, currentMove, nextMove, index) .* (zl .+ prev_zl))


        zl .= preB * b

        #gmres!(solver,opB,b)

        #zl.=solver.x
        zcomp .-= zl

        if norm(zcomp) < error
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


function startup_CPU(PDE_Meta::SchrodingerPDEBase, Grid::SpaceTimeGrid)
    start = zeros(
        Complex{eltype(get_boundary(PDE_Meta, 1))},
        linear_size(get_metadata(Grid)),
        get_total_components(PDE_Meta),
    )
    for (fun, index) in zip(PDE_Meta.ψ₀, 1:get_total_components(PDE_Meta))
        start[:, index] .= evaluate_indexed_elements(Grid, fun)
    end
    start
end

function full_algorithm(Solver::PDESolver, PDE::SchrodingerPDEBase, Grid::SpaceTimeGrid, Ψⁿ)
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
    #plot_inner(currentMove, Grid, 0)
    start_energy = get_diameter(Grid) * sum(abs2.(currentMove), dims = 1)
    for time_step = 1:time_steps

        #Forward
        for component_index = 1:total_components
            update_component!(Solver, PDE, Grid, currentMove, nextMove, component_index)
        end

        #Backward
        for reverse_c_index = total_components:1:-1
            update_component!(Solver, PDE, Grid, currentMove, nextMove, reverse_c_index)
        end


        #display(currentMove)
        if time_step % 100 == 0
            plot_inner(currentMove, Grid, time_step)
            @show get_diameter(Grid) * sum(abs2.(currentMove), dims = 1) - start_energy
            @show time_step
        end
    end
    return nextMove
end

