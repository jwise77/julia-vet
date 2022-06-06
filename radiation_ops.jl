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
Integrate short characteristics for a single 1D cell

    Input
    * I: specific intensity array
    * S: source array
    * chi: opacity array
    * ijk: cell index
    * ray_info: dictionary with number of rays "na", angles (mu = cos phi, theta), and quadature weights

    Output
    * Radiation moments of cell ijk
"""
function integrate_cell_1d(I::AbstractArray, S::AbstractArray, chi::AbstractArray, ijk::Vector{Int}, rays::Dict, omega=1)
    # Calculate ray start/end points
    rstart = ijk' .- rays["dr"]
    rend = ijk' .+ rays["dr"]

### NOTE: just a copy of the 2d (untested) version below. Untouched. Not modified for 1D at all.

    # Arrays for upwind and downwind interpolated quantities
    Ixu = zeros(rays["na"])
    Iyu = zeros(rays["na"])
    Izu = zeros(rays["na"])
    Su = zeros(rays["na"])
    chiu = zeros(rays["na"])
    Sd = zeros(rays["na"])
    chid = zeros(rays["na"])
    
    # Intepolate intensity, opacity, and source functions to those points. Cycling through each Cartesian direction,
    # 1. Create 1D slices of fields for interpolation routine
    # 2. Perform interpolation for all rays
    for dim = 1:2
        # Negative slices
        ijk_m = Dict()
        ijk_zm = Dict()

        for j in 1:2
            ijk_m[j] = (j == dim) ? (ijk[dim] - 1) : Colon()
            ijk_zm[j] = (j == dim) ? (ijk[dim] - 1) : ijk[dim]
        end
        Ixm = view(I, ijk_m[1], ijk_m[2], 1)
        Iym = view(I, ijk_m[1], ijk_m[2], 2)
        Izm = I[ijk_zm[1], ijk_zm[2], 3]
        Sm = view(S, ijk_m[1], ijk_m[2])
        chim = view(chi, ijk_m[1], ijk_m[2])
        # Positive slices
        ijk_p = Dict()
        ijk_zp = Dict()
        for j in 1:2
            ijk_p[j] = (j == dim) ? (ijk[dim] + 1) : Colon()
            ijk_zp[j] = (j == dim) ? (ijk[dim] + 1) : ijk[dim]
        end
        Ixp = view(I, ijk_p[1], ijk_p[2], 1)
        Iyp = view(I, ijk_p[1], ijk_p[2], 2)
        Izp = I[ijk_zp[1], ijk_zp[2], 3]
        Sp = view(S, ijk_p[1], ijk_p[2])
        chip = view(chi, ijk_p[1], ijk_p[2])
        # Execution mask for rays traveling in this direction
        rmask = rays["idir"] .== dim
        neg_rays = rmask .& (rays["isign"] .== -1)
        pos_rays = rmask .& (rays["isign"] .== +1)
        # Interpolation from (-) direction then (+) direction
        otherdim = (dim == 1) ? 2 : 1
        # Upwind
        # Rays traveling in positive direction (starting from -1 plane)
        Ixu[pos_rays] .= bezier_interp_1dv.(Ref(Ixm), rstart[pos_rays,otherdim])
        Iyu[pos_rays] .= bezier_interp_1dv.(Ref(Iym), rstart[pos_rays,otherdim])
        Izu[pos_rays] .= Izm
        Su[pos_rays] .= bezier_interp_1dv.(Ref(Sm), rstart[pos_rays,otherdim])
        chiu[pos_rays] .= bezier_interp_1dv.(Ref(chim), rstart[pos_rays,otherdim])
        # Rays traveling in negative direction (starting from +1 plane)
        Ixu[neg_rays] .= bezier_interp_1dv.(Ref(Ixp), rstart[neg_rays,otherdim])
        Iyu[neg_rays] .= bezier_interp_1dv.(Ref(Iyp), rstart[neg_rays,otherdim])
        Izu[neg_rays] .= Izp
        Su[neg_rays] .= bezier_interp_1dv.(Ref(Sp), rstart[neg_rays,otherdim])
        chiu[neg_rays] .= bezier_interp_1dv.(Ref(chip), rstart[neg_rays,otherdim])
        # Downwind (+1 plane then -1 plane)
        Sd[pos_rays] .= bezier_interp_1dv.(Ref(Sp), rend[pos_rays,otherdim])
        chid[pos_rays] .= bezier_interp_1dv.(Ref(chip), rend[pos_rays,otherdim])
        Sd[neg_rays] .= bezier_interp_1dv.(Ref(Sm), rend[neg_rays,otherdim])
        chid[neg_rays] .= bezier_interp_1dv.(Ref(chim), rend[neg_rays,otherdim])
    end

    # Duplicate source function and opacity at cell center for vectorization
    # Sp and chip are reused from above (which are 2D slices)
    S0 = ones(rays["na"]) * S[ijk[1], ijk[2]]
    chi0 = ones(rays["na"]) * chi[ijk[1], ijk[2]]
    # Calculate interpolated upwind intensity on each ray normal (I cdot n)
    I0 = zeros(rays["na"])
    @. I0 = Ixu * rays["mu"][:,1] + Iyu * rays["mu"][:,2] + Izu * rays["mu"][:,3]
    # Integrate the rays
    Iray = integrate_ray.(I0, Su, S0, Sd, chiu, chi0, chid, rays["ds"], rays["ds"], omega)

    # Calculate moments (J, H, K)
    # Note: Intensity in the Cartesian coordinate directions is the H-vector
    #Inew = [sum(rays["mu"][:,i] .* Iray) for i = 1:3]
    J = sum(rays["w"] .* Iray)
    H = [sum(rays["w"] .* Iray .* rays["mu"][:,i]) for i = 1:3]
    K = reshape([sum(rays["w"] .* Iray .* rays["mu"][:,i] .* rays["mu"][:,j]) 
        for i = 1:3 for j = 1:3]) #, (3,3)
    
    return J, H, K
end

"""
Integrate short characteristics for a single 2D cell

    Input
    * I: specific intensity array
    * S: source array
    * chi: opacity array
    * ijk: 2-vector with cell index
    * ray_info: dictionary with number of rays "na", angles (mu = cos phi, theta), and quadature weights

    Output
    * Radiation moments of cell ijk
"""
function integrate_cell_2d(I::AbstractArray, S::AbstractArray, chi::AbstractArray, ijk::Vector{Int}, rays::Dict, omega=1)
    # Calculate ray start/end points
    rstart = ijk' .- rays["dr"]
    rend = ijk' .+ rays["dr"]

    # Arrays for upwind and downwind interpolated quantities
    Ixu = zeros(rays["na"])
    Iyu = zeros(rays["na"])
    Izu = zeros(rays["na"])
    Su = zeros(rays["na"])
    chiu = zeros(rays["na"])
    Sd = zeros(rays["na"])
    chid = zeros(rays["na"])
    
    # Intepolate intensity, opacity, and source functions to those points. Cycling through each Cartesian direction,
    # 1. Create 1D slices of fields for interpolation routine
    # 2. Perform interpolation for all rays
    for dim = 1:2
        # Negative slices
        ijk_m = Dict()
        ijk_zm = Dict()

        for j in 1:2
            ijk_m[j] = (j == dim) ? (ijk[dim] - 1) : Colon()
            ijk_zm[j] = (j == dim) ? (ijk[dim] - 1) : ijk[dim]
        end
        Ixm = view(I, ijk_m[1], ijk_m[2], 1)
        Iym = view(I, ijk_m[1], ijk_m[2], 2)
        Izm = I[ijk_zm[1], ijk_zm[2], 3]
        Sm = view(S, ijk_m[1], ijk_m[2])
        chim = view(chi, ijk_m[1], ijk_m[2])
        # Positive slices
        ijk_p = Dict()
        ijk_zp = Dict()
        for j in 1:2
            ijk_p[j] = (j == dim) ? (ijk[dim] + 1) : Colon()
            ijk_zp[j] = (j == dim) ? (ijk[dim] + 1) : ijk[dim]
        end
        Ixp = view(I, ijk_p[1], ijk_p[2], 1)
        Iyp = view(I, ijk_p[1], ijk_p[2], 2)
        Izp = I[ijk_zp[1], ijk_zp[2], 3]
        Sp = view(S, ijk_p[1], ijk_p[2])
        chip = view(chi, ijk_p[1], ijk_p[2])
        # Execution mask for rays traveling in this direction
        rmask = rays["idir"] .== dim
        neg_rays = rmask .& (rays["isign"] .== -1)
        pos_rays = rmask .& (rays["isign"] .== +1)
        # Interpolation from (-) direction then (+) direction
        otherdim = (dim == 1) ? 2 : 1
        # Upwind
        # Rays traveling in positive direction (starting from -1 plane)
        Ixu[pos_rays] .= bezier_interp_1dv.(Ref(Ixm), rstart[pos_rays,otherdim])
        Iyu[pos_rays] .= bezier_interp_1dv.(Ref(Iym), rstart[pos_rays,otherdim])
        Izu[pos_rays] .= Izm
        Su[pos_rays] .= bezier_interp_1dv.(Ref(Sm), rstart[pos_rays,otherdim])
        chiu[pos_rays] .= bezier_interp_1dv.(Ref(chim), rstart[pos_rays,otherdim])
        # Rays traveling in negative direction (starting from +1 plane)
        Ixu[neg_rays] .= bezier_interp_1dv.(Ref(Ixp), rstart[neg_rays,otherdim])
        Iyu[neg_rays] .= bezier_interp_1dv.(Ref(Iyp), rstart[neg_rays,otherdim])
        Izu[neg_rays] .= Izp
        Su[neg_rays] .= bezier_interp_1dv.(Ref(Sp), rstart[neg_rays,otherdim])
        chiu[neg_rays] .= bezier_interp_1dv.(Ref(chip), rstart[neg_rays,otherdim])
        # Downwind (+1 plane then -1 plane)
        Sd[pos_rays] .= bezier_interp_1dv.(Ref(Sp), rend[pos_rays,otherdim])
        chid[pos_rays] .= bezier_interp_1dv.(Ref(chip), rend[pos_rays,otherdim])
        Sd[neg_rays] .= bezier_interp_1dv.(Ref(Sm), rend[neg_rays,otherdim])
        chid[neg_rays] .= bezier_interp_1dv.(Ref(chim), rend[neg_rays,otherdim])
    end

    # Duplicate source function and opacity at cell center for vectorization
    # Sp and chip are reused from above (which are 2D slices)
    S0 = ones(rays["na"]) * S[ijk[1], ijk[2]]
    chi0 = ones(rays["na"]) * chi[ijk[1], ijk[2]]
    # Calculate interpolated upwind intensity on each ray normal (I cdot n)
    I0 = zeros(rays["na"])
    @. I0 = Ixu * rays["mu"][:,1] + Iyu * rays["mu"][:,2] + Izu * rays["mu"][:,3]
    # Integrate the rays
    Iray = integrate_ray.(I0, Su, S0, Sd, chiu, chi0, chid, rays["ds"], rays["ds"], omega)

    # Calculate moments (J, H, K)
    # Note: Intensity in the Cartesian coordinate directions is the H-vector
    #Inew = [sum(rays["mu"][:,i] .* Iray) for i = 1:3]
    J = sum(rays["w"] .* Iray)
    H = [sum(rays["w"] .* Iray .* rays["mu"][:,i]) for i = 1:3]
    K = reshape([sum(rays["w"] .* Iray .* rays["mu"][:,i] .* rays["mu"][:,j]) 
        for i = 1:3 for j = 1:3])#, (3,3)
    
    return J, H, K
end

"""
Integrate short characteristics for a single 3D cell

    Input
    * I: specific intensity array
    * S: source array
    * chi: opacity array
    * ijk: 3-vector with cell index
    * ray_info: dictionary with number of rays "na", angles (mu = cos phi, theta), and quadature weights

    Output
    * Radiation moments of cell ijk
"""
function integrate_cell_3d(I::AbstractArray, S::AbstractArray, chi::AbstractArray, ijk::Vector{Int}, rays::Dict, omega=1)
    # Cartestian directions in planes (yz, xz, xy-planes)
    all_idim = ((2,3), (1,3), (1,2))
    # Calculate ray start/end points
    rstart = ijk' .- rays["dr"]
    rend = ijk' .+ rays["dr"]

    # Arrays for upwind and downwind interpolated quantities
    Ixu = zeros(rays["na"])
    Iyu = zeros(rays["na"])
    Izu = zeros(rays["na"])
    Su = zeros(rays["na"])
    chiu = zeros(rays["na"])
    Sd = zeros(rays["na"])
    chid = zeros(rays["na"])

    # Intepolate intensity, opacity, and source functions to those points. Cycling through each Cartesian direction,
    # 1. Create 2D slices of fields for interpolation routine
    # 2. Perform interpolation for all rays
    for dim = 1:3
        # Negative slices
        ijk_m = Dict()
        for j in 1:3
            ijk_m[j] = (j == dim) ? (ijk[dim] - 1) : Colon()
        end
        Ixm = view(I, ijk_m[1], ijk_m[2], ijk_m[3], 1)
        Iym = view(I, ijk_m[1], ijk_m[2], ijk_m[3], 2)
        Izm = view(I, ijk_m[1], ijk_m[2], ijk_m[3], 3)
        Sm = view(S, ijk_m[1], ijk_m[2], ijk_m[3])
        chim = view(chi, ijk_m[1], ijk_m[2], ijk_m[3])
        # Positive slices
        ijk_p = Dict()
        for j in 1:3
            ijk_p[j] = (j == dim) ? (ijk[dim] + 1) : Colon()
        end
        Ixp = view(I, ijk_p[1], ijk_p[2], ijk_p[3], 1)
        Iyp = view(I, ijk_p[1], ijk_p[2], ijk_p[3], 2)
        Izp = view(I, ijk_p[1], ijk_p[2], ijk_p[3], 3)
        Sp = view(S, ijk_p[1], ijk_p[2], ijk_p[3])
        chip = view(chi, ijk_p[1], ijk_p[2], ijk_p[3])
        # Execution mask for rays traveling in this direction
        rmask = rays["idir"] .== dim
        neg_rays = rmask .& (rays["isign"] .== -1)
        pos_rays = rmask .& (rays["isign"] .== +1)
        # Interpolation from -1 planes then +1 planes
        idim = all_idim[dim][1]  # x-direction in plane
        jdim = all_idim[dim][2]  # y-direction
        # Upwind
        # Rays traveling in positive direction (starting from -1 plane)
        Ixu[pos_rays] .= bezier_interp_2d.(Ref(Ixm), rstart[pos_rays,idim], rstart[pos_rays,jdim])
        Iyu[pos_rays] .= bezier_interp_2d.(Ref(Iym), rstart[pos_rays,idim], rstart[pos_rays,jdim])
        Izu[pos_rays] .= bezier_interp_2d.(Ref(Izm), rstart[pos_rays,idim], rstart[pos_rays,jdim])
        Su[pos_rays] .= bezier_interp_2d.(Ref(Sm), rstart[pos_rays,idim], rstart[pos_rays,jdim])
        chiu[pos_rays] .= bezier_interp_2d.(Ref(chim), rstart[pos_rays,idim], rstart[pos_rays,jdim])
        # Rays traveling in negative direction (starting from +1 plane)
        Ixu[neg_rays] .= bezier_interp_2d.(Ref(Ixp), rstart[neg_rays,idim], rstart[neg_rays,jdim])
        Iyu[neg_rays] .= bezier_interp_2d.(Ref(Iyp), rstart[neg_rays,idim], rstart[neg_rays,jdim])
        Izu[neg_rays] .= bezier_interp_2d.(Ref(Izp), rstart[neg_rays,idim], rstart[neg_rays,jdim])
        Su[neg_rays] .= bezier_interp_2d.(Ref(Sp), rstart[neg_rays,idim], rstart[neg_rays,jdim])
        chiu[neg_rays] .= bezier_interp_2d.(Ref(chip), rstart[neg_rays,idim], rstart[neg_rays,jdim])
        # Downwind (+1 plane then -1 plane)
        Sd[pos_rays] .= bezier_interp_2d.(Ref(Sp), rend[pos_rays,idim], rstart[pos_rays,jdim])
        chid[pos_rays] .= bezier_interp_2d.(Ref(chip), rend[pos_rays,idim], rstart[pos_rays,jdim])
        Sd[neg_rays] .= bezier_interp_2d.(Ref(Sm), rend[neg_rays,idim], rstart[neg_rays,jdim])
        chid[neg_rays] .= bezier_interp_2d.(Ref(chim), rend[neg_rays,idim], rstart[neg_rays,jdim])
    end

    # Duplicate source function and opacity at cell center for vectorization
    # Sp and chip are reused from above (which are 2D slices)
    S0 = ones(rays["na"]) * S[ijk[1], ijk[2], ijk[3]]
    chi0 = ones(rays["na"]) * chi[ijk[1], ijk[2], ijk[3]]
    # Calculate interpolated upwind intensity on each ray normal (I cdot n)
    I0 = zeros(rays["na"])
    @. I0 = Ixu * rays["mu"][:,1] + Iyu * rays["mu"][:,2] + Izu * rays["mu"][:,3]
    # Integrate the rays
    Iray = integrate_ray.(I0, Su, S0, Sd, chiu, chi0, chid, rays["ds"], rays["ds"], omega)

    # Calculate moments (J, H, K)
    # Note: Intensity in the Cartesian coordinate directions is the H-vector
    #Inew = [sum(rays["mu"][:,i] .* Iray) for i = 1:3]
    J = sum(rays["w"] .* Iray)
    H = [sum(rays["w"] .* Iray .* rays["mu"][:,i]) for i = 1:3]
    K = reshape([sum(rays["w"] .* Iray .* rays["mu"][:,i] .* rays["mu"][:,j]) 
        for i = 1:3 for j = 1:3], (3,3))
    
    return J, H, K
end
