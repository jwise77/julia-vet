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
    * (u)pwind, (d)ownwind, (p)oint
* S: Source function dictionary, same abbrevations
* ds: path length on the (u)pwind / (d)ownwind rays
* omega: Bezier shape parameter (0.5 = parabolic, 1 = linear)
Output: Intensity
"""
function integrate_ray(I0::Number, Su::Number, Sp::Number, Sd::Number, 
    chiu::Number, chip::Number, chid::Number,
    ds_u::Number, ds_d::Number, omega=1)

    # Equations 13-14: Optical depth between upwind/downwind points and control point
    # Calculate control point and then optical depths
    chic_u = bezier_controlm(chiu, chip, chid, omega)
    chic_d = bezier_control(chiu, chip, chid, omega)
    dtau_u = (ds_u / 3.0) * (chiu + chic_u + chip)
    dtau_d = (ds_d / 3.0) * (chid + chic_d + chip)

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
    result = a * Su + b * Sp + c * Sd + d * I0
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
function integrate_cell(I::AbstractArray, S::AbstractArray, chi::AbstractArray, ijk::Vector{Int}, rays::Dict, omega=1)
    # Calculate ray start/end points
    rstart = ijk' .- rays["dr"]
    rend = ijk' .+ rays["dr"]

    # Arrays for upwind and downwind interpolated quantities
    Iu = zeros(rays["na"])
    Su = zeros(rays["na"])
    chiu = zeros(rays["na"])
    Id = zeros(rays["na"])
    Sd = zeros(rays["na"])
    chid = zeros(rays["na"])

    # Intepolate intensity, opacity, and source functions to those points. Cycling through each Cartesian direction,
    # 1. Create 2D slices of fields for interpolation routine
    # 2. Perform interpolation for all rays
    for dim = 1:3
        # Negative slices
        ijk_m = (:,:,:)
        ijk_m[dim] = ijk_m[dim] - 1
        Ixm = view(I, ijk_m[1], ijk_m[2], ijk_m[3], 1)
        Iym = view(I, ijk_m[1], ijk_m[2], ijk_m[3], 2)
        Izm = view(I, ijk_m[1], ijk_m[2], ijk_m[3], 3)
        Sm = view(S, ijk_m[1], ijk_m[2], ijk_m[3])
        chim = view(chi, ijk_m[1], ijk_m[2], ijk_m[3])
        # Positive slices
        ijk_p = (:,:,:)
        ijk_p[dim] = ijk_p[dim] + 1
        Ixp = view(I, ijk_p[1], ijk_p[2], ijk_p[3], 1)
        Iyp = view(I, ijk_p[1], ijk_p[2], ijk_p[3], 2)
        Izp = view(I, ijk_p[1], ijk_p[2], ijk_p[3], 3)
        Sp = view(S, ijk_p[1], ijk_p[2], ijk_p[3])
        chip = view(chi, ijk_p[1], ijk_p[2], ijk_p[3])
        # Execution mask for rays traveling in this direction
        rmask = rays["idir"] == dim
        neg_rays = rmask .& (rays["isign"] == -1)
        pos_rays = rmask .& (rays["isign"] == +1)
        # Interpolation (upwind) from -1 planes then +1 planes
        Ixu[neg_rays] = bezier_interp_2d.(Ixm, rstart[neg_rays])
        Iyu[neg_rays] = bezier_interp_2d.(Iym, rstart[neg_rays])
        Izu[neg_rays] = bezier_interp_2d.(Izm, rstart[neg_rays])
        Su[neg_rays] = bezier_interp_2d.(Sm, rstart[neg_rays])
        chiu[neg_rays] = bezier_interp_2d.(chim, rstart[neg_rays])
        Ixu[pos_rays] = bezier_interp_2d.(Ixp, rstart[pos_rays])
        Iyu[pos_rays] = bezier_interp_2d.(Iyp, rstart[pos_rays])
        Izu[pos_rays] = bezier_interp_2d.(Izp, rstart[pos_rays])
        Su[pos_rays] = bezier_interp_2d.(Sp, rstart[pos_rays])
        chiu[pos_rays] = bezier_interp_2d.(chip, rstart[pos_rays])
        # Downwind
        Sd[neg_rays] = bezier_interp_2d.(Sm, rend[neg_rays])
        chid[neg_rays] = bezier_interp_2d.(chim, rend[neg_rays])
        Sd[pos_rays] = bezier_interp_2d.(Sp, rend[pos_rays])
        chid[pos_rays] = bezier_interp_2d.(chip, rend[pos_rays])
    end

    # Duplicate source function and opacity at cell center for vectorization
    Sp = ones(rays["da"]) * S[ijk[1], ijk[2], ijk[3]]
    chip = ones(rays["da"]) * chi[ijk[1], ijk[2], ijk[3]]
    # Calculate interpolated upwind intensity on each ray normal (I cdot n)
    I0 = Ixu .* rays["mu"][:,1] .+ Iyu .* rays["mu"][:,2] .+ Izu .* rays["mu"][:,3]
    # Integrate the rays
    Iray = integrate_ray.(I0, Su, Sp, Sd, chiu, chip, chid, rays["ds"], rays["ds"], omega)

    # Update intensity and moments (J, H, K)
    Inew = [sum(rays["mu"][:,i] .* Iray) for i = 1:3]
    J = sum(rays["w"] .* Iray)
    H = [sum(rays["w"] .* Iray .* rays["mu"][:,i]) for i = 1:3]
    K = reshape([sum(rays["w"] .* Iray .* rays["mu"][:,i] .* rays["mu"][:,j]) 
        for i = 1:3 for j = 1:3], (3,3))
    
    return Inew, J, H, K
end