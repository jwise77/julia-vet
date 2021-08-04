import Distributions: Uniform
include("math_routines.jl")

N = 5

# Zeros
# one = Array{Float64}(undef, N)
# two = Array{Float64}(undef, N, N)
# three = Array{Float64}(undef, N, N, N)

# Random
rfn = Uniform()
one = rand(rfn, N)
two = rand(rfn, N,N)
three = rand(rfn, N,N,N)

fp = grid_deriv(one)
fp2 = grid_deriv(two)
fp3 = grid_deriv(three)

#
# Test interpolation
#
# 1D interpolation points (for 2D RT problems along adjacent edges)
nx = 3  # number of points
xx = rand(Uniform(1,N),nx)
fx = interp_cubic(one, fp, xx)
println(fx)

# 2D interpolation points (for 3D RT problems on adjacent faces)
xx2 = rand(Uniform(1,N),nx,2)

