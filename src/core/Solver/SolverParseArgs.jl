"""
    This function is used for argument parsing and construction of default structures. This function don't relies on a specific backend.
    
    Return:
               τ times coefficients   unique values of prev   values           unique σ values    
    N-dim Grid,TimeCollection vector,TimeMultipliers vector,Time Steps count, Set of σ values
"""
function parse_input_parameters(::Type{T},PDE::SchrodingerPDE,
    time_order::Int = 2,
    time_composition_substeps::Int = 1,
    time_composition_index::Int = 1,kwargs...)::Tuple{SpaceTimeGrid,Vector{T},Vector{T},Int,Set}

    #Dimensionality
    expected_dims = ndims(PDE)

    #Validations and startup calculations...
    h, N, dim_ranges, τ, time_steps = parse_input_dim_args(T,
        expected_dims,
        PDE;
        kwargs...)

    #Abstract Mesh with data...
    Grid = CreateGrid(dim_ranges..., τ, h..., N...)

    #Initialization of Time Discretization

    TimeDiscretization = get_time_discretization(T,
        time_order,
        time_composition_substeps,
        time_composition_index)

    #Initialization of Space Discretization

    TimeMultipliers = τ * get_coefficients(TimeDiscretization)
    TimeCollection = τ * collect(TimeDiscretization)
    σset = Set(get_σ(PDE))


    Grid,TimeCollection,TimeMultipliers,time_steps,σset
end