include("types.jl")

function set_compute_backend(new_compute_backend::String)
    if !(new_compute_backend in ("OpenCL", "CUDA", "CPU","CPUMultiThread"))
        throw(ArgumentError("Invalid backend: \"$(new_compute_backend)\""))
    end

    # Set it in our runtime values, as well as saving it to disk
    @set_preferences!("compute_backend" => new_compute_backend)
    @info("New backend set; restart your Julia session for this change to take effect!")
end

const compute_backend = @load_preference("compute_backend","CPU")

@show compute_backend