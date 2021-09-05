using Debugger
import Distributions: Uniform
include("angular_gauss_quad.jl")
include("bezier-h20.jl")
include("radiation_ops.jl")

struct grid
    rank::Int
    dims::AbstractArray
    S::AbstractArray
    I::AbstractArray
    B::AbstractArray
    chi::AbstractArray
    rho::AbstractArray
    T::AbstractArray
    J::AbstractArray
    H::AbstractArray
    K::AbstractArray
end

function test_bezier()
    
    N = 6
    
    # Zeros
    # one = Array{Float64}(undef, N)
    # two = Array{Float64}(undef, N, N)
    # three = Array{Float64}(undef, N, N, N)
    
    # Random
    rfn = Uniform()
    one = rand(rfn, N)
    two = rand(rfn, N,N)
    three = rand(rfn, N,N,N)
    
    #fp = grid_deriv(one)
    #fp2 = grid_deriv(two)
    #fp3 = grid_deriv(three)
    
    #
    # Test interpolation
    #
    # 1D interpolation points (for 2D RT problems along adjacent edges)
    nx = 3  # number of points
    xa = rand(Uniform(2,N-1),nx)
    x = xa[1]
    fx = bezier_interp_1dv(one, x)
    println(fx)
    
    # 2D interpolation points (for 3D RT problems on adjacent faces)
    xx2 = rand(Uniform(2,N-1),nx,2)
    xx = xx2[:,1]
    yy = xx2[:,2]
    xi2 = floor.(Int, xx2)
    xi = xi2[:,1]
    yi = xi2[:,2]
    fxy = bezier_interp_2d.(two, xx, yy)
    println(fxy)
end

function initialize_regular_grid(rank::Int, N::Int)
    dims = N * ones(Int, rank)
    empty_grid = zeros(Tuple(dims))
    empty_vec_grid = zeros(Tuple([dims;rank]))
    empty_tensor_grid = zeros(Tuple([dims;rank;rank]))
    result = grid(rank, dims,
        copy(empty_grid), # S
        copy(empty_vec_grid), # I
        copy(empty_grid), # B
        copy(empty_grid), # chi
        copy(empty_grid), # rho
        copy(empty_grid), # (T)emperature
        copy(empty_grid), # J
        copy(empty_vec_grid), # H
        copy(empty_tensor_grid) # K
        )
    return result
end

function test_ray()
    I0 = 0.0
    S = Dict("u"=>1.0, "p"=>0.1, "d"=>0.1)
    chi = Dict("u"=>0.1, "p"=>0.5, "d"=>1.0)
    Inew = integrate_ray(I0, S["u"], S["p"], S["d"], chi["u"], chi["p"], chi["d"], 1, 1, 1)
end

function test_cell()
    rank = 3
    N = 10
    nmu = 4
    omega = 0.5  # 0.5 (parabolic) -> 1.0 (linear) interpolation
    grid = initialize_regular_grid(rank, N)
    grid.chi .= 0.1
    grid.S[1,:,:] .= 1.0
    rays = calculate_ray_info(nmu)
    ijk = [2,2,2]
    J, H, K = integrate_cell(grid.I, grid.S, grid.chi, ijk, rays, omega)
end
