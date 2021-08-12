using Debugger
import Distributions: Uniform
#include("general_interpolation.jl")
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
#fx = bezier_interp_1dv(one, x)
#println(fx)

# 2D interpolation points (for 3D RT problems on adjacent faces)
xx2 = rand(Uniform(2,N-1),nx,2)
xx = xx2[:,1]
yy = xx2[:,2]
xi2 = floor.(Int, xx2)
xi = xi2[:,1]
yi = xi2[:,2]