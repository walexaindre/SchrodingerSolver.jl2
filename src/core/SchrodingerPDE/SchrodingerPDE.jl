include("SchrodingerPDETypes.jl")


@inline function get_space_dimension(PDE_Meta::SchrodingerPDE)
    size(PDE_Meta.boundaries, 2)
end

@inline function get_total_components(PDE_Meta::SchrodingerPDE)
    size(PDE_Meta.σ, 1)
end

@inline function eval_component(PDE_Meta::SchrodingerPDENonPolynomic, index::Int, arguments...)
    PDE_Meta.f[index](arguments...)
end

@inline function eval_initial_condition(PDE_Meta::SchrodingerPDE, index::Int, arguments...)
    PDE_Meta.ψ₀[index](arguments...)
end

@inline function eval_optimized_potential(PDE_Meta::SchrodingerPDEPolynomic,prevIteration,nextIteration,index)
    PDE_Meta.N(prevIteration,nextIteration,index)
end

@inline function eval_field(PDE_Meta::SchrodingerPDE, arguments...)
    PDE_Meta.F(arguments...)
end

@inline function get_σ(PDE_Meta::SchrodingerPDE, index::Int)
    PDE_Meta.σ[:,index]
end

@inline function get_boundary(PDE_Meta::SchrodingerPDE, index::Int)
    PDE_Meta.boundaries[:, index]
end

@inline function get_end_time(PDE_Meta::SchrodingerPDE)
    PDE_Meta.T
end