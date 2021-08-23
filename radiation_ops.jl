include("angular_gauss_quad.jl")
# Methods to operator on radiation in grids and rays

"""
Return the opacity for a given state
"""
function compute_opacity(rho::Number)
    # For testing purposes, set them equal
    chi = rho
    return chi
end

"""
Blackbody radiation
B(nu,T) = \frac{2h\nu^3}{c^2} \frac {1}{e^{h\nu / kT}-1}
For testing, set pre-factor to one.
Input x = h\nu/kt
"""
function Bnu(x::Number)
    return 1.0 / (exp(x) - 1.0)
end 
function Bnu(T::Number, nu::Number)
    x = 4.79924e-11 * T / nu  # h\nu / kT
    return 1.0 / (exp(x) - 1.0)
end

"""
Update source function
    For now, do nothing
"""
function source_update!(source)
    return nothing
end

"""
Integrate RTE along short characteristic (Hennenicker+ 2020, Equation 12)
* I_ijk = a_ijk * S_upwind + b_ijk * S_point + c_ijk * S_downwind + d_ijk * I_upwind

Input
* I0: original intensity
* chi: opacity dictionary
    * (u)pwind, (d)ownwind, (c)ontrol_(u)pwind, (c)ontrol_(d)ownwind, (p)oint
* S: Source function dictionary, same abbrevations
* ds: path length between between (u)pwind / (d)ownwind points and control point
* omega: Bezier shape parameter (0.5 = parabolic, 1 = linear)
Output: Intensity
"""
function integrate_ray(I0::Number, S::Dict, chi::Dict, ds_u::Number, ds_d::Number, omega=1)
    # Equations 13-14: Optical depth between upwind/downwind points and control point
    dtau_u = (ds_u / 3.0) * (chi["u"] + chi["c_u"] + chi["p"])
    dtau_d = (ds_d / 3.0) * (chi["d"] + chi["c_d"] + chi["p"])

    # Integrals of the RTE w.r.t. optical depth (Equation 12+)
    e0 = 1 - exp(-dtau_u)
    e1 = dtau_u - e0
    e2 = dtau_u^2 - 2 * e1

    # Coefficients (Equation 12+)
    a = e0 + (omega - 2) * e1 / dtau_u + (1 - omega) * e2 / dtau_u^2
    b = ((1 - omega) * dtau_u + (2 - omega) * dtau_d) * e1 / (dtau_u * dtau_d) +
        (omega - 1) * (dtau_u + dtau_d) * e2 / (dtau_u^2 * dtau_d)
    c = (omega - 1) * e1 / dtau_d + (1 - omega) * e2 / (dtau_u * dtau_d)
    d = exp(-dtau_u)

    # Integrate
    result = a * S["u"] + b * S["p"] + c * S["d"] + d * I0
    return result
end

"""
Integrate short characteristics for a single cell

    Input
    * I: specific intensity array
    * S: source array
    * chi: opacity array
    * ijk: 3-vector with cell index (ones in places where the rank is higher than the dimension)
    * ray_info: dictionary with number of rays "na", angles (mu = cos phi, theta), and quadature weights

    Output
    * New intensity and associated moments of cell ijk
"""
function integrate_cell(I::AbstractArray, S::AbstractArray, chi::AbstractArray, ijk::Vector{Int}, rays::Dict)
    # Calculate ray start/end points
    rstart = ijk' .- rays["dr"]
    rend = ijk' .+ rays["dr"]

    # Intepolate intensity, opacity, and source functions to those points. Cycling through each Cartesian direction,
    # 1. Create 2D slices of fields for interpolation routine
    # 2. Perform interpolation for all rays
    for dim = 1:3
        # Negative slices
        ijk_u = (:,:,:)
        ijk_u[dim] = ijk_u[dim] - 1
        Iu = view(I, ijk_u[1], ijk_u[2], ijk_u[3])
        Su = view(S, ijk_u[1], ijk_u[2], ijk_u[3])
        chiu = view(chi, ijk_u[1], ijk_u[2], ijk_u[3])
        # Positive slices
        ijk_d = (:,:,:)
        ijk_d[dim] = ijk_d[dim] + 1
        Id = view(I, ijk_d[1], ijk_d[2], ijk_d[3])
        Sd = view(S, ijk_d[1], ijk_d[2], ijk_d[3])
        chid = view(chi, ijk_d[1], ijk_d[2], ijk_d[3])
        # Execution mask for rays traveling in this direction
        rmask = rays["idir"] == dim
        rdir = sign.(rays["dr"][:,rays["idir"]])  # pointing in -1 or +1 direction
        # Interpolation
        #Iu_rays = bezier_interp_2d.(Id, rstart)
    end

    # Update intensity and moments (J, H, K)

end