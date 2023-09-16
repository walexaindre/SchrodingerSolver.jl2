export generate_kernel_B, generate_kernel_C,generate_preconditioner_B,assembly_metakernel,get_time_steps

function count_nonzero_elements(Mesh::MetaMesh1D)
    3 * linear_size(Mesh)
end

function count_nonzero_elements(Mesh::MetaMesh2D)
    5 * linear_size(Mesh)
end

function count_nonzero_elements(Mesh::MetaMesh3D)
    7 * linear_size(Mesh)
end


function generate_kernel_B(::Type{T}, Mesh::MetaMesh1D,τ::T,σ::T,hx::T) where {T<:AbstractFloat}
    nonzero = count_nonzero_elements(Mesh)
    n = linear_size(Mesh)

    sqhx = (τ * σ) / (hx .^ 2) 

    diag = (4im|>Complex{T})-(2|>T)*(sqhx)

    I = zeros(Int64, nonzero)
    J = zeros(Int64, nonzero)
    V = zeros(Complex{T}, nonzero)

    @threads for index = 1:n
        ix = linear_indexing(index,Mesh)
        ax,bx=step_x(ix,Mesh)

        ix=linear_indexing(ix,Mesh)
        ax=linear_indexing(ax,Mesh)
        bx=linear_indexing(bx,Mesh)

        #@show [ixiy;axiy;bxiy;ixay;ixby]
        @views I[3*(index-1)+1:3*(index)] .= index
        @views J[3*(index-1)+1:3*(index)] .= [ix;ax;bx]
        @views V[3*(index-1)+1:3*(index)] .= [diag;sqhx;sqhx]

    end

#@show J

    #S[I[k], J[k]] = V[k].
    sparse(I, J, V, n, n)
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
        @views I[5*(index-1)+1:5*(index)] .= index
        @views J[5*(index-1)+1:5*(index)] .= [ixiy;axiy;bxiy;ixay;ixby]
        @views V[5*(index-1)+1:5*(index)] .= [diag;sqhx;sqhx;sqhy;sqhy]

    end

#@show J

    #S[I[k], J[k]] = V[k].
    sparse(I, J, V, n, n)
end

function generate_kernel_B(::Type{T}, Mesh::MetaMesh3D,τ::T,σ::T,hx::T,hy::T,hz::T) where {T<:AbstractFloat}
    nonzero = count_nonzero_elements(Mesh)
    n = linear_size(Mesh)

    sqhx = (τ * σ) / (hx .^ 2) 
    sqhy = (τ * σ) / (hy .^ 2)
    sqhz = (τ * σ) / (hz .^ 2)

    diag = (4im|>Complex{T})-(2|>T)*(sqhx+sqhy+sqhz)

    I = zeros(Int64, nonzero)
    J = zeros(Int64, nonzero)
    V = zeros(Complex{T}, nonzero)

    @threads for index = 1:n
        ix,iy,iz = linear_indexing(index,Mesh)
        ay,by=step_x(iy,Mesh)
        ax,bx=step_y(ix,Mesh)
        az,bz=step_z(iz,Mesh)

        ixiyiz=linear_indexing(ix,iy,iz,Mesh)
        axiyiz=linear_indexing(ax,iy,iz,Mesh)
        bxiyiz=linear_indexing(bx,iy,iz,Mesh)
        ixayiz=linear_indexing(ix,ay,iz,Mesh)
        ixbyiz=linear_indexing(ix,by,iz,Mesh)
        ixiyaz=linear_indexing(ix,iy,az,Mesh)
        ixiybz=linear_indexing(ix,iy,bz,Mesh)
        

        #@show [ixiy;axiy;bxiy;ixay;ixby]
        @views I[7*(index-1)+1:7*(index)] .= index
        @views J[7*(index-1)+1:7*(index)] .= [ixiyiz;axiyiz;bxiyiz;ixayiz;ixbyiz;ixiyaz;ixiybz]
        @views V[7*(index-1)+1:7*(index)] .= [diag;sqhx;sqhx;sqhy;sqhy;sqhz;sqhz]

    end

#@show J

    #S[I[k], J[k]] = V[k].
    sparse(I, J, V, n, n)
end


function generate_kernel_C(::Type{T}, Mesh::MetaMesh1D,τ::T,σ::T,hx::T) where {T<:AbstractFloat}
    return generate_kernel_B(T,Mesh,-τ,σ,hx)
end

function generate_kernel_C(::Type{T}, Mesh::MetaMesh2D,τ::T,σ::T,hx::T,hy::T) where {T<:AbstractFloat}
    return generate_kernel_B(T,Mesh,-τ,σ,hx,hy)
end

function generate_kernel_C(::Type{T}, Mesh::MetaMesh3D,τ::T,σ::T,hx::T,hy::T,hz::T) where {T<:AbstractFloat}
    return generate_kernel_B(T,Mesh,-τ,σ,hx,hy,hz)
end

function generate_preconditioner_B(A::SparseMatrixCSC)
    return ilu(A)
end


function assembly_metakernel(::Type{CPUBackend},Mesh::MetaMesh,B::SparseMatrixCSC,C::SparseMatrixCSC,iluB::IncompleteLU.ILUFactorization)
    opB = LinearOperator(B)
    opC = LinearOperator(C)

    n = linear_size(Mesh)

    #opB = LinearOperator(Complex{DataType}, n, n, true, true, b_call)
    preB = LinearOperator(eltype(B),n,n,false,false,(y,x)-> IncompleteLU.backward_substitution!(iluB,IncompleteLU.forward_substitution!(y,iluB,x)))
    
    MetaKernel(opB,opC,preB)
end

function get_time_steps(PDE::SchrodingerPDEBase, Grid::SpaceTimeGrid)
    PDE.T / Grid.τ
end