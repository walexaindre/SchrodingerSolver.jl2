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

        #Used stencil positions at stencil by droptolerance.
        stencil_positions = findall(in(int_index),stencil)

        used_stencil_idx = stencil[stencil_positions]


        nnz_stencil_positions = length(stencil_positions)

        I = Array{Int64}(undef,rows*nnz_stencil_positions)

        J = Array{Int64}(undef,rows*nnz_stencil_positions)

        V = repeat(col[used_stencil_idx]|>Array,cols)


        #@show "Stencil used"
        #display(stencil_used_idx)
        @show "Stencil pos"
        display(stencil_positions)
        @show "int ndex"
        display(int_index)

        @threads for idx in 1:length(Mesh)
            I[(nnz_stencil_positions * (idx - 1) + 1):(nnz_stencil_positions * idx)] .= get_linear_stencil(idx,
            1,
            stencil_max_depth,
            Mesh)[stencil_positions]
    
            #Column is always the same.
            J[(nnz_stencil_positions * (idx - 1) + 1):(nnz_stencil_positions * idx)] .= idx
        end
    return I,J,V
end
