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

function bezier_interp_1d(f, x, omega=1, coeff_only=false)
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
    omegai = 1 ./ (1 .- (f[idxp] - f[idx0]) / (f[idx0] - f[idxm]))
    omegap = 1 .+ 1 ./ (1 .- (f[idx0] - f[idxp]) / (f[idxp] - f[idx0]))
    omega0 = 0.5 * (omegai + omegap)
    # Restrict omega to be within the range [omegai, omegap].
    # If outside, set to be the arthimetic mean of omegai and omegap.
    om = omega .* ones(Float64, size(x))
    mask = (omega .< omegai) .| (omega .> omegap)
    om[mask] = omega0[mask]

    # Interpolation coefficients. (B.14-16) setting dx = 1
    atilde = (om - 1) .* t .- (om - 1) .* t.^2
    btilde = 1 .- (2*om - 1) .* t .+ 2*(om - 1) .* t.^2
    ctilde = om .* t + (1 - om) .* t.^2

    # Evaluate the Bezier polynomial with the coefficients
    if coeff_only
        return atilde, btilde, ctilde
    else
        result = atilde .* f[idxm] .+ btilde .* f[idx0] .+ ctilde .* f[idx1]
        return result
    end
end