export generate_kernel_B

function count_nonzero_elements(Mesh::MetaMesh1D)
    3 * linear_size(Mesh)
end

function count_nonzero_elements(Mesh::MetaMesh2D)
    5 * linear_size(Mesh)
end

function count_nonzero_elements(Mesh::MetaMesh3D)
    7 * linear_size(Mesh)
end

function generate_kernel_B(::Type{T}, Mesh::MetaMesh2D,τ::T,σ::T,hx::T,hy::T) where {T<:AbstractFloat}
    nonzero = count_nonzero_elements(Mesh)
    n = linear_size(Mesh)

    sqhx = (τ * σ) / (hx .^ 2) 
    sqhy = (τ * σ) / (hy .^ 2)

    diag = (4im|>Complex{T})-(2|>T)*(sqhx+sqhy)

    I = zeros(Int64, nonzero)
    J = zeros(Int64, nonzero)
    V = zeros(Complex{T}, nonzero)

    @threads for index = 1:n
        ix,iy = linear_indexing(index,Mesh)
        ay,by=step_x(iy,Mesh)
        ax,bx=step_y(ix,Mesh)

        ixiy=linear_indexing(ix,iy,Mesh)
        axiy=linear_indexing(ax,iy,Mesh)
        bxiy=linear_indexing(bx,iy,Mesh)
        ixay=linear_indexing(ix,ay,Mesh)
        ixby=linear_indexing(ix,by,Mesh)

        #@show [ixiy;axiy;bxiy;ixay;ixby]
        @views I[5*(index-1)+1:5*(index)] .= [index;index;index;index;index]
        @views J[5*(index-1)+1:5*(index)] .= [ixiy;axiy;bxiy;ixay;ixby]
        @views V[5*(index-1)+1:5*(index)] .= [diag;sqhx;sqhx;sqhy;sqhy]

    end

#@show J

    #S[I[k], J[k]] = V[k].
    sparse(I, J, V, n, n)
end