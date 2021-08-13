#=
Grid interpolation routines for upstream / downstream ray values
=#

struct field_deriv
    ndims::Int
    size::Tuple
    left::Array{Array}
    right::Array{Array}
end

"""
Returns the left- and right-handed differences of a field
"""
function grid_deriv(f, dx=1)
    nd = ndims(f)
    ns = size(f)
    #vec_size = Tuple(Base.Iterators.flatten((nd,ns)))
    fprime = field_deriv(nd, Tuple(ns), [zeros(ns) for i=1:nd], [zeros(ns) for i=1:nd])
    for i in 1:nd
        shift = zeros(Int,nd)
        shift[i] = 1
        #fprime.left[:,:,:,i]  .= f .- circshift(f,shift)
        #fprime.right[:,:,:,i] .= circshift(f,-shift) .- f
        fprime.left[i] .= (f .- circshift(f,shift)) ./ dx
        fprime.right[i] .= (circshift(f,-shift) .- f) ./ dx
    end
    return fprime
end

# Interpolation functions in 1D
# Perform multi-D interpolation in sweeps

"""
Nearest neighbor interpolation
"""
    function interp_zero(f,fp,x,dx=1)
    idx = round.(Int, x/dx)
    return f[idx]
end

"""
1D Linear Interpolation
"""
    function interp_linear(f,fp,x,dx=1)
    t = mod.(x,dx)
    idx = floor.(Int, x/dx)
    result = f[idx] .+ t.*(f[idx.+1] .- f[idx])/dx
    return result
end

"""
Hayek+ 2010, Equation B.3
Harmonic mean of left- and right-handed differences to remove wiggles and overshoots in strong gradients
"""
    function fp_harm(fp,dim=1)
    fpnew = zeros(size(fp.left[dim]))
    pp = fp.left[dim] .* fp.right[dim] .> 0
    fpnew[pp] = (fp.left[dim][pp] .* fp.right[dim][pp]) ./ (0.5 * (fp.left[dim][pp] .+ fp.right[dim][pp]))
    return fpnew
end

"""
1D Monotonic quadratic (Hayek+ 2010, Appendix B, equation B.5)
"""
function interp_quad(f,fp,x,dx=1,dim=1)
    nx = size(x)
    t = mod.(x,dx)
    idx = floor.(Int, x/dx)

    left = t .< 0.5
    right = t .> 0.5
    alpha = zeros(nx)
    beta = zeros(nx)
    gamma = zeros(nx)
    delta = zeros(nx)
    # Interpolation on the left-hand side
    alpha[left] = 1 .- t[left].^2
    beta[left] = t[left].^2
    gamma[left] = t[left] .* (1 .- t[left]) * dx
    # Interpolation on the right-hand side
    alpha[right] = (1 .- t[right]).^2
    beta[right] = t[right] .* (2 .- t[right])
    delta[right] = t[right] .* (t[right] .- 1) * dx
    fpnew = fp_harm(fp,dim)
    result = alpha .* f[idx] .+ beta .* f[idx.+1] .+ gamma .* fpnew[idx] .+ delta .* fpnew[idx.+1]
    return result
end

"""
1D Monotonic cubic (Hayek+ 2010, Appendix B)
"""
function interp_cubic(f,fp,x,dx=1,dim=1)
    nx = size(x)
    t = mod.(x, dx)
    idx = floor.(Int, x/dx)
    alpha = 2*t.^3 .- 3*t.^2 .+ 1
    beta = 3*t.^2 .- 2*t.^3
    gamma = (t.^3 .- 2*t.^2 .+ t)*dx
    delta = (t.^3 .- t.^2)*dx
    fpnew = fp_harm(fp,dim)
    result = alpha .* f[idx] .+ beta .* f[idx.+1] .+ gamma .* fpnew[idx] .+ delta .* fpnew[idx.+1]
    return result
end

"""
General interpolation routine: Only 1D currently
Interpolate positions x (rank order) from field f (rank order)
Assume that x are in units of cell widths (dx = 1)
"""
    function interpolate(f, x, order=1)
    if order == 0
        interp_fn = interp_zero
    elseif order == 1
        interp_fn = interp_linear
    elseif order == 2
        interp_fn = interp_quad
    elseif order == 3
        interp_fn = interp_cubic
    end
    # Left and right derivatives
    fp = grid_deriv(f)
    if ndims(f) > 1
        error("Multi-dimensional interpolation not implemented.")
    end
    result = interp_fn(f, fp, x)
    return result
end
