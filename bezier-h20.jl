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

function bezier_interp_1dv(f, x, omega=1, return_coeff=false)
    # Version that takes entire array (row or column)
    # coeff_only
    # * true = returns interpolation coefficients only. Useful for 2D.
    # * false = return interpolated value
    idx0 = floor.(Int, x)
    idxm = idx0 .- 1
    idxp = idx0 .+ 1
    xi = Float64.(idx0)
    # Interpolation point (B.11)
    t = x .- xi
    # Control point (B.12)
    fc = f[idx0] .+ 0.5 * (omega .* (f[idxp] .- f[idx0]) .+ 
        (1 .- omega) .* (f[idx0] - f[idxm]))

    # Shape parameters to maintain monotonicity (B.17-18)
    # Note: Not sure if Julia will deal with divide by zeroes nicely
    omegai = 1 ./ (1 .- (f[idxp] .- f[idx0]) ./ (f[idx0] .- f[idxm]))
    omegap = 1 .+ 1 ./ (1 .- (f[idx0] .- f[idxm]) ./ (f[idxp] .- f[idx0]))
    omega_min = min.(omegai, omegap)
    omega_max = max.(omegai, omegap)
    omega0 = 0.5 * (omega_min .+ omega_max)
    # Restrict omega to be within the range [omegai, omegap].
    # If outside, set to be the arthimetic mean of omegai and omegap.
    om = omega .* ones(Float64, size(x))
    mask = (omega .< omega_min) .| (omega .> omega_max)
    #om[mask] = omega0[mask]
    om[mask] = max.(min.(omega_min[mask], 0.5), 1.0)

    # Interpolation coefficients. (B.14-16) setting dx = 1
    atilde = (om .- 1) .* t .- (om .- 1) .* t.^2
    btilde = 1 .- (2 .* om .- 1) .* t .+ 2 .* (om .- 1) .* t.^2
    ctilde = om .* t .+ (1 .- om) .* t.^2

    # Evaluate the Bezier polynomial with the coefficients
    result = atilde .* f[idxm] .+ btilde .* f[idx0] .+ ctilde .* f[idxp]
    if return_coeff
        coeff = Dict("a"=>atilde, "b"=>btilde, "c"=>ctilde)
        result = (result, coeff)
    end
    return result
end

function bezier_interp_1d(fm, f0, fp, x, omega=1, return_coeff=false)
    # Version that takes three points (j-1, j, j+1) and an x-coordinate for
    # interpolation instead of an array.
    # Same operations as the full array version.
    t = mod.(x,1)
    fc = f0 .+ 0.5 * (omega .* (fp .- f0) .+ (1 .- omega) .* (f0 .- fm))

    omegai = 1 ./ (1 .- (fp .- f0) ./ (f0 .- fm))
    omegap = 1 .+ 1 ./ (1 .- (f0 .- fm) ./ (fp .- f0))
    omega_min = min.(omegai, omegap)
    omega_max = max.(omegai, omegap)
    omega0 = 0.5 * (omegai .+ omegap)

    om = omega .* ones(Float64, size(x))
    mask = (omega .< omega_min) .| (omega .> omega_max)
    #om[mask] = omega0[mask]
    om[mask] = max.(min.(omega_min[mask], 0.5), 1.0)

    atilde = (om .- 1) .* t .- (om .- 1) .* t.^2
    btilde = 1 .- (2 .* om .- 1) .* t .+ 2 .* (om .- 1) .* t.^2
    ctilde = om .* t + (1 .- om) .* t.^2

    result = atilde .* fm .+ btilde .* f0 .+ ctilde .* fp
    if return_coeff
        coeff = Dict("a"=>atilde, "b"=>btilde, "c"=>ctilde)
        result = (result, coeff)
    end
    return result
end

function get_indices(f, ix, iy)
    # Does what I expect f[ix,iy] would do
    # Returns an array of size (ix) instead of (ix,iy)
    return [f[ix[i],iy[i]] for i in eachindex(ix,iy)]
end

function bezier_interp_2d(f, x, y, omega=1)
    if ndims(f) != 2
        error("Field f must be 2D.")
    end

    # Central and adjacent indices
    idx0 = floor.(Int, x)
    idxm = idx0 .- 1
    idxp = idx0 .+ 1
    idy0 = floor.(Int, y)
    idym = idy0 .- 1
    idyp = idy0 .+ 1

    # Interpolated values (fy) and coefficients (cy) for j-1, j, and j+1 rows
    fxm, cxm = bezier_interp_1dv(f[:,idym], x, omega, true)
    fx0, cx0 = bezier_interp_1dv(f[:,idy0], x, omega, true)
    fxp, cxp = bezier_interp_1dv(f[:,idyp], x, omega, true)

    # Interpolation values (finterp) and coefficients (cx) along the y-axis 
    # at the x-value
    fy, cy = bezier_interp_1d(fxm, fx0, fxp, y, omega, true)

    # 2D interpolation with all coefficients (Equation C.1)
    result = 
        cy["a"] .* cxm["a"] .* get_indices(f, idxm, idym) + 
        cy["a"] .* cxm["b"] .* get_indices(f, idx0, idym) +
        cy["a"] .* cxm["c"] .* get_indices(f, idxp, idym) +
        cy["b"] .* cx0["a"] .* get_indices(f, idxm, idy0) + 
        cy["b"] .* cx0["b"] .* get_indices(f, idx0, idy0) +
        cy["b"] .* cx0["c"] .* get_indices(f, idxp, idy0) +
        cy["c"] .* cxp["a"] .* get_indices(f, idxm, idyp) + 
        cy["c"] .* cxp["b"] .* get_indices(f, idx0, idyp) +
        cy["c"] .* cxp["c"] .* get_indices(f, idxp, idyp)
    
    return result
end