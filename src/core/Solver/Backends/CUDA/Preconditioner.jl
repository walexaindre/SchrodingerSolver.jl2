export drop
function drop(M::CuSparseMatrixCSR, Mesh::MetaMesh, τ = 0.0001::AbstractFloat, stencil_max_depth::Int=3, r_tol::AbstractFloat = 700 * eps(real(eltype(M))),
        a_tol::AbstractFloat = 700 * eps(real(eltype(M)))) 
    
        rows,cols=size(M)

        e₁ = CUDA.zeros(eltype(M),rows)
        e₁[1:1] .= 1.0
        col, _ = gmres(M,e₁,atol=a_tol,rtol=r_tol)
        #Get the number of nonzero entries in the inverse of the first column.
        nnz_entries = count(x->!isapprox(abs2(x),0,atol=a_tol),col)

        #Prune the inverse with the desired drop tolerance τ.
        bitindex = abs2.(col) .> τ
        #Get indexes of non pruned elements.
        int_index = findall(bitindex)|>Array

        #Stencil with desired depth
        stencil = get_linear_stencil(1,1,stencil_max_depth,Mesh)
        
        #Valid stencil positions at int_index.
        stencil_used_idx = findall(in(stencil),int_index)
        
        #Used stencil positions at stencil.
        stencil_positions = findall(in(int_index),stencil)

        nnz_stencil_positions = length(stencil_positions)

        I = Array{Int64}(undef,rows*nnz_stencil_positions)

        J = Array{Int64}(undef,rows*nnz_stencil_positions)

        V = repeat(col[int_index][stencil_positions]|>Array,cols)

        display(col[int_index][stencil_positions]|>Array)
        display(col[int_index]|>Array)
        display(col[1:13]|>Array)
        #gmres

        display(get_linear_stencil(1,
        1,
        stencil_max_depth,
        Mesh)[stencil_positions])

        display(get_linear_stencil(2,
        1,
        stencil_max_depth,
        Mesh)[stencil_positions])

        @threads for idx in 1:length(Mesh)
            I[(nnz_stencil_positions * (idx - 1) + 1):(nnz_stencil_positions * idx)] .= get_linear_stencil(idx,
            1,
            stencil_max_depth,
            Mesh)[stencil_positions]
    
            #Column is always the same.
            J[(nnz_stencil_positions * (idx - 1) + 1):(nnz_stencil_positions * idx)] .= idx
        end


        #Iterator lenght step. (Columns are independent. Construction of preconditioner can be done in parallel)
        itlen = length(stencil_used_idx)

    ##@warn
    return I,J,V
end
