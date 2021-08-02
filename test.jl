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
nx = 3
xx = rand(Uniform(1,N),nx)
