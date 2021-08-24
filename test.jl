using Debugger
include("angular_gauss_quad.jl")
include("radiation_ops.jl")

function test_bezier()
    import Distributions: Uniform
    include("bezier-h20.jl")
    
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
    fxy = bezier_interp_2d(two, xx, yy)
    println(fxy)
end

function initialize_regular_grid(rank::Int, N::Int)
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
    dims = N * ones(rank, Int)
    empty_grid = zeros(dims)
    empty_vec_grid = zeros([dims;rank])
    empty_tensor_grid = zeros([dims;rank;rank])
    result = grid(rank, dims,
        empty_grid, # S
        empty_vec_grid, # I
        empty_grid, # B
        empty_grid, # chi
        empty_grid, # rho
        empty_grid, # (T)emperature
        empty_grid, # J
        empty_vec_grid, # H
        empty_tensor_grid # K
        )
    return result
end

function test_ray()
    I0 = 0.0
    S = Dict("u"=>1.0, "p"=>0.1, "d"=>0.0)
    chi = Dict("u"=>0.1, "p"=>0.5, "d"=>1.0)
    Inew = integrate_ray(I0, S["u"], S["p"], S["d"], chi["u"], chi["p"], chi["d"], 1, 1, 1)
end

function test_cell()
    rank = 3
    N = 10
    nmu = 4
    omega = 1  # 0.5 (parabolic) -> 1.0 (linear) interpolation
    grid = initialize_regular_grid(3, N)
    rays = calculate_ray_info(nmu)
    ijk = [2,2,2]
    Inew, J, H, K = integrate_cell(grid.I, grid.S, grid.chi, ijk, rays, omega)
end
