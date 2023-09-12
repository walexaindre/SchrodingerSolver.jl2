include("types.jl")


function get_space_dimension(PDE_Meta::SchrodingerPDEBase)
    size(PDE_Meta.boundaries, 2)
end

function get_total_components(PDE_Meta::SchrodingerPDEBase)
    size(PDE_Meta.σ, 1)
end

function eval_component(PDE_Meta::SchrodingerPDENonPolynomic, index::Int, arguments...)
    PDE_Meta.f[index](arguments...)
end

function eval_initial_condition(PDE_Meta::SchrodingerPDEBase, index::Int, arguments...)
    PDE_Meta.ψ₀[index](arguments...)
end

function eval_optimized_potential(PDE_Meta::SchrodingerPDEPolynomic,prevIteration,nextIteration)
    PDE_Meta.N(prevIteration,nextIteration)
end

function eval_field(PDE_Meta::SchrodingerPDEBase, arguments...)
    PDE_Meta.F(arguments...)
end

function get_σ(PDE_Meta::SchrodingerPDEBase, index::Int)
    PDE_Meta.σ[index]
end

function get_boundary(PDE_Meta::SchrodingerPDEBase, index::Int)
    PDE_Meta.boundaries[:, index]
end

function get_end_time(PDE_Meta::SchrodingerPDEBase)
    PDE_Meta.T
end

function compute()

   @static if (compute_backend=="CUDA")
    println("Computing with CUDA")
   else
    println("Computing with CPU")
   end 
    
end

