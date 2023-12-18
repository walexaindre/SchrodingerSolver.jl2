"Obtain the count of non-zero entries for square matrices using specified diagonals. The main diagonal is assigned the index 0. '-idx' denotes diagonals below the main diagonal, while 'idx' refers to diagonals above the main diagonal."
function nnz_diag_entries(dim::Int,used_diags::Array{Int,1})
    length(used_diags)*dim-sum(abs,used_diags)
end
