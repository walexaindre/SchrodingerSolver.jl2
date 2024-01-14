export drop

function drop(M::CuSparseMatrixCSR, Mesh::MetaMesh, τ = 0.0001::AbstractFloat;
        stencil_min_depth::Int = 1, stencil_max_depth::Int = 3,
        r_tol::AbstractFloat = 700 * eps(real(eltype(M))),
        a_tol::AbstractFloat = 700 * eps(real(eltype(M))))
    rows, _ = size(M)

    e₁ = CUDA.zeros(eltype(M), rows)
    e₁[1:1] .= 1.0

    col, _ = gmres(M, e₁, atol = a_tol, rtol = r_tol) #Ax=e₁ (First column inverse)

    #Prune the inverse with the desired drop tolerance τ.
    bit_index = abs2.(col) .> τ

    #Get indexes of non pruned elements.
    int_index = findall(bit_index) |> Array

    #Stencil with desired depth
    full_stencil = get_linear_stencil(1, stencil_min_depth, stencil_max_depth, Mesh)

    #Undropped stencil positions.
    undropped_stencil_positions = findall(in(int_index), full_stencil)

    #Here the not eliminated entries are chosen.
    used_stencil_idx = full_stencil[undropped_stencil_positions]

    nnz_stencil_positions = length(undropped_stencil_positions)

    #Checking a_ij==a_ji
    if 1:nnz_stencil_positions != undropped_stencil_positions
        @warn "Problem with matrix entries (a_ij!=a_ji). Use a different drop tolerance..."
    end

    #e₁ values at used stencil positions.
    ê₁ = col[used_stencil_idx] |> Array

    #Row index depends on prunned stencil for every column.
    I = Array{Int64}(undef, rows * nnz_stencil_positions)

    #Column is always the same even nnz_stencil_positions.
    J = Array{Int64}(undef, rows * nnz_stencil_positions)

    #Values are always the same for every column.
    V = Array{Int64}(undef, rows * nnz_stencil_positions)

    #IJV Construction can be done in parallel.
    @threads for idx in 1:length(Mesh)
        b = nnz_stencil_positions * idx
        a = b - nnz_stencil_positions + 1

        pos_to_modify = a:b

        #Rows are generated by prunned stencil expansion at current position. 
        I[pos_to_modify] .= get_linear_stencil(idx,
            stencil_min_depth,
            stencil_max_depth,
            Mesh)[undropped_stencil_positions]

        #Column is always the same.
        J[pos_to_modify] .= idx

        #Values are always the same.
        V[pos_to_modify] .= ê₁
    end
    return I, J, V
end