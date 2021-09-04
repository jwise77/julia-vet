# 1D and 2D Bezier interpolation from Hennicker et al. (2020)
# Gives a better and more thorough description than the other papers.
#
# We only need the 1D and 2D cases for short characteristics because we
# need the quantities on the ray intersections with the cell edges 
# (2D domain) or cell faces (3D domain). No need for linear interpolation
# routines because of the control parameter omega. In Appendix B of 
# Hennicker+ (2020), we will only use the [x_i, x_(i+1)] case because we'll 
# floor the values.
#
# Also see bezier-hayek.jl
#
# =========================================================================
#
# Input:
# * f: 1D or 2D field
# * x: actual x-positions in units of cell widths
# * y: (optional) same but for y-positions
# * omega: Bezier control parameter between linear (=1) and quadratic (=0.5)
#
# Assumption: all cell widths are equal
#
# Output: interpolated quantity

"""
Restrict the interpolation shape parameter to preserve monotonicity
"""
function restrict_omega(fm::Number, f0::Number, fp::Number, omega::Number, plus_side::Bool)
    dfp = fp - f0
    dfm = f0 - fm
    # Prevent divide by zero
    dfm = max(abs(dfm), 1e-20)
    dfp = max(abs(dfp), 1e-20)
    one_dfp = 1 - dfp / dfm
    one_dfm = 1 - dfm / dfp
    # Prevent divide by zero
    one_dfm = max(abs(one_dfm), 1e-20)
    one_dfp = max(abs(one_dfp), 1e-20)
    if plus_side
        omegai = 1 / one_dfp
        omegap = 1 + 1 / one_dfm
    else
        omegai = 1 + 1 / one_dfp
        omegap = 1 / one_dfm
    end
    # Restrict omega to range of [omegai,omegap] and [0.5,1]
    omega_min = max(min(omegai, omegap), 0.5)
    omega_max = min(max(omegai, omegap), 1.0)
    om = min(max(omega, omega_min), omega_max)
    return om
end

"""
Hennenick+ (2020), Appendix B.2: Interval [x_i, x_{i+1}]
* For quantity q interpolations, we always assume that x is in this range.
* Version that takes entire array (row or column)
* coeff_only
    * true = returns interpolation coefficients only. Useful for 2D.
    * false = return interpolated value
"""
function bezier_interp_1dv(f::AbstractArray, x::Number, omega=1, return_coeff=false)
    if ndims(f) != 1
        error("f must be one-dimensional")
    end
    idx0 = floor(Int, x)
    idxm = max(idx0-1, 1)
    idxp = min(idx0+1, size(f)[1])
    omega = (idx0 == idxm || idx0 == idxp) ? 1.0 : omega  # linear if boundary cell

    # Interpolation point (B.11)
    xi = Float64(idx0)
    t = x - xi

    # Shape parameters to maintain monotonicity (B.17-18)
    om = restrict_omega(f[idxm], f[idx0], f[idxp], omega, true)

    # Interpolation coefficients. (B.14-16) setting dx = 1
    atilde = (om - 1) * t - (om - 1) * t^2
    btilde = 1 - (2 * om - 1) * t + 2 * (om - 1) * t^2
    ctilde = om * t + (1 - om) * t^2

    # Control point (B.12)
    fc = f[idx0] + 0.5 * (om * (f[idxp] - f[idx0]) + 
        (1 - om) * (f[idx0] - f[idxm]))

    # Evaluate the Bezier polynomial with the coefficients
    result = atilde * f[idxm] + btilde * f[idx0] + ctilde * f[idxp]
    if return_coeff
        coeff = Dict("a"=>atilde, "b"=>btilde, "c"=>ctilde, "fc"=>fc)
        result = (result, coeff)
    end
    return result
end

"""
Hennenick+ (2020), Appendix B.2: Interval [x_i, x_{i+1}]
* For quantity q interpolations, we always assume that x is in this range.
* Version that takes three points (j-1, j, j+1) and an x-coordinate for interpolation instead of an array.
* Same operations as the full array version.
"""
function bezier_interp_1d(fm::Number, f0::Number, fp::Number, x::Number, omega=1, return_coeff=false)
    t = mod(x,1)

    om = restrict_omega(fm, f0, fp, omega, true)
    dfp = fp - f0
    atilde = (om - 1) * t - (om - 1) * t^2
    btilde = 1 - (2 * om - 1) * t + 2 * (om - 1) * t^2
    ctilde = om * t + (1 - om) * t^2
    fc = f0 + 0.5 * (om * (fp - f0) + (1 - om) * (f0 - fm))

    result = atilde * fm + btilde * f0 + ctilde * fp
    if return_coeff
        coeff = Dict("a"=>atilde, "b"=>btilde, "c"=>ctilde, "fc"=>fc)
        result = (result, coeff)
    end
    return result
end

"""
Hennenick+ (2020), Appendix B.1: Interval [x_{i-1}, x_i]
* Only used for interpolation along a ray in the upstream portion, in particular the opacity and associated control point.
* Version that takes three points (j-1, j, j+1) and an x-coordinate for interpolation instead of an array.
* Same operations as the full array version.
"""
function bezier_interp_1dm(fm::Number, f0::Number, fp::Number, x::Number, omega=1, return_coeff=false)
    t = mod(x,1)

    om = restrict_omega(fm, f0, fp, omega, false)

    atilde = 1 + (om - 2) * t + (1 - om) * t^2
    btilde = (3 - 2 * om) * t + 2 * (om - 1) * t^2
    ctilde = (om - 1) * t - (om - 1) * t^2
    fc = f0 + 0.5 * (om * (f0 - fm) + (1 - om) * (fp - f0))

    result = atilde * fm + btilde * f0 + ctilde * fp
    if return_coeff
        coeff = Dict("a"=>atilde, "b"=>btilde, "c"=>ctilde, "fc"=>fc)
        result = (result, coeff)
    end
    return result
end

"""
Hennenick+ (2020), Appendix B.1 (Equations B.5-6) [i-1,i] interval
* Return function at control point. Depends on function values for a uniform grid.
"""
function bezier_controlm(fm::Number, f0::Number, fp::Number, omega=1)
    om = restrict_omega(fm, f0, fp, omega, false)
    fc = f0 + 0.5 * (om * (fp - f0) + (1 - om) * (f0 - fm))
    return fc
end

"""
Hennenick+ (2020), Appendix B.2 (Equations B.17-18) [i,i+1] interval
* Return function at control point. Depends on function values for a uniform grid.
"""
function bezier_control(fm::Number, f0::Number, fp::Number, omega=1)
    om = restrict_omega(fm, f0, fp, omega, true)
    fc = f0 + 0.5 * (om * (fp - f0) + (1 - om) * (f0 - fm))
    return fc
end

"""
Does what I expect f[ix,iy] would do
Returns an array of size (ix) instead of (ix,iy)
"""
function get_indices(f::AbstractMatrix, ix::AbstractArray, iy::AbstractArray)
    return [f[ix[i],iy[i]] for i in eachindex(ix,iy)]
end

function bezier_interp_2d(f::AbstractArray, x::Number, y::Number, omega=1)
    if ndims(f) != 2
        error("Field f must be 2D.")
    end
    if (x < 1) || (x > size(f)[1])
        error("Coordinate x ($x) out of bounds")
    end
    if (y < 1) || (y > size(f)[2])
        error("Coordinate y ($y) out of bounds")
    end

    # Central and adjacent indices
    idx0 = floor(Int, x)
    idy0 = floor(Int, y)
    if idx0 > 1
        idxm = idx0 - 1
        omegam = omega
    else
        idxm = 1
        omegam = 1.0  # linear interpolation
    end
    if idx0 < size(f)[1]
        idxp = idx0 + 1
        omegap = omega
    else
        idxp = size(f)[1]
        omegap = 1.0  # linear interpolation
    end

    idym = max(idy0-1, 1)
    idyp = min(idy0+1, size(f)[2])

    # Interpolated values (fy) and coefficients (cy) for j-1, j, and j+1 rows
    fxm, cxm = bezier_interp_1dv.(Ref(f[:,idym]), x, omegam, true)
    fx0, cx0 = bezier_interp_1dv.(Ref(f[:,idy0]), x, omega, true)
    fxp, cxp = bezier_interp_1dv.(Ref(f[:,idyp]), x, omegap, true)

    # Interpolation values (finterp) and coefficients (cx) along the y-axis 
    # at the x-value
    omegay = (idy0 == idym || idy0 == idyp) ? 1.0 : omega
    fy, cy = bezier_interp_1d.(fxm, fx0, fxp, y, omegay, true)

    # 2D interpolation with all coefficients (Equation C.1)
    result = 
        cy["a"] * cxm["a"] * f[idxm, idym] + 
        cy["a"] * cxm["b"] * f[idx0, idym] +
        cy["a"] * cxm["c"] * f[idxp, idym] +
        cy["b"] * cx0["a"] * f[idxm, idy0] + 
        cy["b"] * cx0["b"] * f[idx0, idy0] +
        cy["b"] * cx0["c"] * f[idxp, idy0] +
        cy["c"] * cxp["a"] * f[idxm, idyp] + 
        cy["c"] * cxp["b"] * f[idx0, idyp] +
        cy["c"] * cxp["c"] * f[idxp, idyp]
    
    return result
end
