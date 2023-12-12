include("SpaceDiscretizationTypes.jl")
export get_space_discretization_coefficients
export get_sparsematrix_A, get_sparsematrix_D

function get_space_discretization_coefficients(::Type{T},
    order::Int)::SecondDerivativeCoefficients{T} where {T <: AbstractFloat}
    if order == 2
        α = 0 // 1
        β = 0 // 1

        a = 1 // 1
        b = 0 // 1
        c = 0 // 1

    elseif order == 4
        α = 1 // 10
        β = 0 // 1

        a = (4 // 3) * (1 - α)
        b = (1 // 3) * (-1 + 10 * α)
        c = 0 // 1
    elseif order == 6
        α = 2 // 11
        β = 0 // 1

        a = 12 // 11
        b = 3 // 11
        c = 0 // 1
    elseif order == 8
        α = 344 // 1179
        β = (38 * α - 9) // 214

        a = (696 - 1191 * α) // 428
        b = (2454 * α - 294) // 535
        c = (1179 * α - 344) // 2140

    elseif order == 10
        α = 334 // 899
        β = 43 // 1798

        a = 1065 // 1798
        b = 1038 // 899
        c = 79 // 1798
    else
        throw(BoundsError([2, 4, 6, 8, 10], order))
    end

    SecondDerivativeCoefficients{T}(T(a), T(b), T(c), T(α), T(β))
end

function get_A_format_IJV(::Type{T},
    Mesh::MetaMesh,
    order::Int) where {T <: AbstractFloat}
    space_discretization = get_space_discretization_coefficients(T, order)

    α = space_discretization.α
    β = space_discretization.β
    count = 0
    value = [T(ndims(Mesh))]
    start_depth = 0
    end_depth = 0

    if (α != 0)
        count += 2
        start_depth = 1
        end_depth = 1
        append!(value,fill(α, 2 * ndims(Mesh)))
    end
    if (β != 0)
        count += 2
        end_depth = 2
        append!(value,fill(β, 2 * ndims(Mesh)))
    end

    if (α == 0 && β != 0)
        start_depth = 2
    end

    count = count * ndims(Mesh) + 1

    space_usage = count * length(Mesh)

    I = zeros(Int64, space_usage) #row idx
    J = zeros(Int64, space_usage) #column idx
    V = repeat(value, length(Mesh))

    @threads for idx in 1:length(Mesh)
        I[(count * (idx - 1) + 1):(count * idx)] .= idx
        if (start_depth != 0)
            J[(count * (idx - 1) + 1):(count * idx)] .= get_linear_stencil(idx,
                start_depth,
                end_depth,
                Mesh)
        else
            J[(count * (idx - 1) + 1):(count * idx)] .= idx
        end

        #V[(count * (idx - 1) + 1):(count * idx)] .= value
    end
    I, J, V
end

@inline function get_sparsematrix_A(::Type{T},
    Mesh::MetaMesh1D,
    order::Int) where {T <: AbstractFloat}
    sparse(get_A_format_IJV(T, Mesh, order)...)
end

@inline function get_A_csr(::Type{T},
    Mesh::MetaMesh1D,
    order::Int) where {T <: AbstractFloat}
    sparsecsr(get_A_format_IJV(T, Mesh, order)...)
end

@inline function get_sparsematrix_A(::Type{T},
    Mesh::MetaMesh2D,
    order::Int) where {T <: AbstractFloat}

    MeshX = MetaMesh1D(Mesh.N)
    MeshY = MetaMesh1D(Mesh.M)

    Ax = get_sparsematrix_A(T, MeshX, order)

    Ay = get_sparsematrix_A(T, MeshY, order)

    kron(Ay, Ax)
end

function Δ_format_IJV(n::Int, coefficient::T, h::T, lsize::Int) where {T <: AbstractFloat}
    if !(0 <= n <= 2)
        throw(BoundsError([0, 1, 2], n))
    end

    if n == 0
        hmul = 1
    elseif n == 1
        hmul = 4
    else
        hmul = 9
    end

    diag = (-2 * coefficient) / (hmul * h .^ 2)
    nodiag = (coefficient) / (hmul * h .^ 2)

    values = [diag, nodiag, nodiag]

    mesh = MetaMesh1D(lsize)

    space_usage = 3 * length(mesh)

    I = zeros(Int64, space_usage) #row idx
    J = zeros(Int64, space_usage) #column idx
    V = repeat(values, length(mesh))

    @threads for idx in 1:length(mesh)
        section = (3 * (idx - 1) + 1):(3 * idx)
        I[section] .= idx
        J[section] .= get_linear_stencil(idx, n + 1, n + 1, mesh)
    end
    I, J, V
end

function Δ(n::Int, coefficient::T, h::T, lsize::Int) where {T <: AbstractFloat}
    sparse(Δ_format_IJV(n, coefficient, h, lsize)...)
end

function Δ_csr(n::Int, coefficient::T, h::T, lsize::Int) where {T <: AbstractFloat}
    sparsecsr(Δ_format_IJV(n, coefficient, h, lsize)...)
end

function get_sparsematrix_D(Grid::SpaceTimeGrid1D, order::Int)
    coeff = get_space_discretization_coefficients(typeof(Grid.hx), order)

    a = coeff.a
    b = coeff.b
    c = coeff.c

    x_size = size(Grid)
    hx = get_space_steps(Grid)

    Δ(0, a, hx, x_size) + Δ(1, b, hx, x_size) + Δ(2, c, hx, x_size)
end

function get_sparsematrix_D(Grid::SpaceTimeGrid2D{T}, order::Int) where {T <: AbstractFloat}
    coeff = get_space_discretization_coefficients(typeof(Grid.hx), order)

    a = coeff.a
    b = coeff.b
    c = coeff.c

    x_size, y_size = size(Grid)
    hx, hy = get_space_steps(Grid)

    Ix = sparse(T(1.0)I, x_size, x_size)
    Iy = sparse(T(1.0)I, y_size, y_size)

    Dx = Δ(0, a, hx, x_size) + Δ(1, b, hx, x_size) + Δ(2, c, hx, x_size)
    Dy = Δ(0, a, hy, y_size) + Δ(1, b, hy, y_size) + Δ(2, c, hy, y_size)
    
    kron(Iy, Dx) + kron(Dy, Ix)
end

function get_sparsematrix_D(Grid::SpaceTimeGrid3D{T}, order::Int) where {T <: AbstractFloat}
    coeff = get_space_discretization_coefficients(typeof(Grid.hx), order)

    a = coeff.a
    b = coeff.b
    c = coeff.c

    x_size, y_size, z_size = size(Grid)
    hx, hy, hz = get_space_steps(Grid)

    Ix = sparse(T(1.0)I, x_size, x_size)
    Iy = sparse(T(1.0)I, y_size, y_size)
    Iz = sparse(T(1.0)I, z_size, z_size)

    Dx = Δ(0, a, hx, x_size) + Δ(1, b, hx, x_size) + Δ(2, c, hx, x_size)
    Dy = Δ(0, a, hy, y_size) + Δ(1, b, hy, y_size) + Δ(2, c, hy, y_size)
    Dz = Δ(0, a, hz, z_size) + Δ(1, b, hz, z_size) + Δ(2, c, hz, z_size)

    kron(Iz, Iy, Dx) + kron(Iz, Dy, Ix) + kron(Dz, Iy, Ix)
end

function Δₕ(Grid::SpaceTimeGrid3D{T}, order::Int) where {T <: AbstractFloat}
end