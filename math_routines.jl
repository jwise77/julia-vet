# Calculate Gaussian angular quadrature: method of discrete ordinates
# (Bruls et al. 1999, Appendix B)
#
# mu = cos(theta)
#
# nmu=2: point for each quadrant.
# Athena-RT (Davis et al. 2012) uses up to nmu=12
# Results in a na = nmu*(nmu+2) total number of rays

using Printf, StructArrays

function ang_weights_mu(nmu, verbose=false)

    if nmu % 2 == 1
        println("Number of polar angles must be even.  Setting n_mu = $nmu -> $(nmu+1)")
        nmu += 1
    elseif nmu == 1
        # Trivial answer for one point
        mu = [sqrt(1/3)]
        w = [1.0]
        theta = [acos(mu[1])]
        return w, mu, theta
    end

    na = nmu * (nmu+2)

    if verbose
        println("We'll calculate $nmu polar angles, resulting in $na total rays.")
    end

    # Calculate angles and weights on north hemisphere only because
    # they're symmetric

    mu1_sq = 1/(3*(nmu-1))

    # Use Newton-Raphson (secant) method to find W1
    l = collect(1:nmu/2-1)
    delta = 2/(nmu-1)

    # Define function and its derivative
    f(w1) = sum(sqrt.(w1^2 .+ (l.-1).*delta)) .- (nmu-2)/3
    fp(w1) = sum(w1 ./ sqrt.(w1^2 .+ (l.-1).*delta))

    tol = 1e-4
    maxiter = 5
    last = 1.0 # first estimate
    next = 4/(3*(nmu-1))  # second estimate ( always will be too large
                          # ). Taken from Appendix B in Bruls et
                          # al. (1999).

    count = 0
    while abs(next - last) > tol && count < maxiter
        last = next
        next = last - f(last)/fp(last)
        #@printf "new guess[%d] = %g -> %g; delta = %g\n" count last next abs(last-next)
        count += 1
    end
    W1_sq = next^2
    ########################################################################
    # END Newton-Raphson
    ########################################################################
    
    dmu = (2 / (nmu-2)) * (1 - 3*mu1_sq)
    dW = 2/(nmu-1)
    j = collect(1:nmu/2)
    mu_sq = mu1_sq .+ (j.-1)*dmu
    W_sq = W1_sq .+ (j.-1)*dW

    # Calculate mu, theta, and w for northern hemisphere

    bigW = sqrt.(W_sq)
    wN = zeros(div(nmu,2))
    wN[1] = bigW[1]
    if nmu > 2
        wN[2:end-1] = bigW[2:end-1] - bigW[1:end-2]
        wN[end] = 1.0 - sum(wN[1:end-1])
    end
    muN = sqrt.(mu_sq)
    thetaN = acos.(muN)

    # Apply equatorial symmetry. Points lie on circles with a co-latitude theta
    w = zeros(nmu)
    mu = zeros(nmu)
    theta = zeros(nmu)
    
    mid = Int(nmu/2)
    # South
    theta[1:mid] = -thetaN
    mu[1:mid] = muN
    w[1:mid] = wN
    
    # North
    theta[mid+1:end] = reverse(thetaN)
    mu[mid+1:end] = reverse(muN)
    w[mid+1:end] = reverse(wN)

    if verbose
        println("mu = $mu")
        println("theta = $theta")
        println("w = $w")
        @printf "sum(w) = %f\n" sum(w)
    end

    return w, mu, theta
    
end

function bezier_control(s0, sm, sp, dxm, dxp)
    # Calculate control point in a Bezier curve
    dx_sum = dxm + dxp
    s0prime = (dxp / dx_sum) * (s0 - sm) / dxm + (dxm / dx_sum) * (sp - s0) / dxp
    sc = s0 - 0.5 * dxm * s0prime
    return sc
end

function bezier_interp_coef(s0, sm, sp, dtaum, dtaup, linear=false)
    # Kunasz & Auer (1988), Equations (7a) - (9c)
    # Hayek et al. (2010), Appendix A, psim = Psi_d to conform with Kunasz88 notation
    u0(x) = 1 - exp(-x)
    u1(x) = x - u0(x)
    u2(x) = x^2 - 2*u1(x)
    if linear
        psi0 = u1(dtaum) / dtaum
        psim = u0(dtaum) - psi0
        psip = 0.0
    else
        # second-order interpolation
        psi0 = ((dtaum + dtaup) * u1(dtaum) - u2(dtaum)) / (dtaum * dtaup)
        psim = u0(dtaum) + (u2(dtaum) - (dtaup + 2*dtaum) * u1(dtaum)) / (dtaum * (dtaum + dtaup))
        psip = (u2(dtaum) - dtaum * u1(dtaum)) / (dtaup * (dtaum + dtaup))
        # Hayek et al. (2010) Appendix A
        # Limit overshooting
        sc = bezier_control(s0, sm, sp, dtaum, dtaup)
        if max(sm,s0,sp) == s0
            # if S0 (cell center) is the extremum, choose S0 as the control point in the Bezier curve
            psi0 = 2 * dtaum * u1(dtaum) - u2(dtaum) / dtaum^2
            psim = u0(dtaum) - psi0
            psip = 0.0
        else if sc < min(sm,s0) || sc > max(sm,s0)
            # if Sc is outside the data range (i.e. overshooting), choose the upwind point as the control point
            psi0 = u2(dtaum) / dtaum^2
            psim = u0(dtaum) - psi0
            psip = 0.0
        end
    end
    return psi0, psim, psip
end

function bezier_interp_tau(chi0, chim, chip, dsm, dsp)
    # Bezier interpolation for optical depth (Equation A.3 in Hayek+ 2010)
    chic = bezier_control(chi0, chim, chip, dsm, dsp)
    # Limit overshooting
    if max(chim, chi0, chip) == chi0
        chic = chi0
    else if chic < min(chim, chi0) || chic > max(chim, chi0)
        chic = chim
    end
    ds = dsm + dsp
    # Interpolation
    dtau = ds * (chim + chi0 + chic) / 3
    return dtau
end

#=
Grid interpolation routines for upstream / downstream ray values
=#

struct field_deriv{T<:Real}
    ndims::Int
    size::Tuple
    left::Array{T}
    right::Array{T}
end

function grid_deriv(f, dx=1)
    nd = ndims(f)
    ns = size(f)
    vec_size = Tuple(Base.Iterators.flatten((ns,nd)))
    fprime = field_deriv(nd, ns, zeros(vec_size), zeros(vec_size))
    for i in 1:nd
        shift = zeros(Int,nd)
        shift[i] = 1
        fprime.left[:,:,i]  = f - circshift(f,shift)
        fprime.right[:,:,i] = circshift(f,-shift) - f
    end
    fprime.left /= dx
    fprime.right /= dx
    return fprime
end

# Interpolation functions in 1D
# Perform multi-D interpolation in sweeps

function interp_zero(f,fp,x,dx=1)
    # Nearest neighbor
    idx = floor.(Int, x/dx)
    return f[idx]
end

function interp_linear(f,fp,x,dx=1)
    # 1D Linear Interpolation (sweeps)
    t = mod.(x,dx)
    idx = floor.(Int, x/dx)
    result = f[idx] .+ t.*(f[idx+1] .- f[idx-1])/dx
    return result
end

function fp_harm(fp)
    # Hayek+ 2010, Equation B.3
    # Harmonic mean of left- and right-handed differences to remove wiggles and overshoots in strong gradients
    fpnew = zeros(fp.size)
    pp = fp.left * fp.right .> 0
    fpnew[pp] = (fp.left[pp] .* fp.right[pp]) ./ (0.5 * (fp.left[pp] .+ fp.right[pp]))
    return fpnew
end

function interp_quad(f,fp,x,dx=1)
    # 1D Monotonic quadratic (Hayek+ 2010, Appendix B, equation B.5)
    nx = size(x)
    t = mod.(x,dx)
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
    fpnew = fp_harm(fp)
    result = alpha .* f[idx] + beta .* f[idx+1] + gamma .* fpnew[idx] + delta .* fpnew[idx+1]
    return result
end

function interp_cubic(f,fp,x,dx=1)
    # 1D Monotonic cubic (Hayek+ 2010, Appendix B)
    nx = size(x)
    t = mod.(x, dx)
    idx = floor.(Int, x/dx)
    alpha = 2*t.^3 .- 3*t.^2 .+ 1
    beta = 3*t.^2 .- 2*t.^3
    gamma = (t.^3 .- 2*t.^2 .+ t)*dx
    delta = (t.^3 .- t.^2)*dx
    fpnew = fp_harm(fp)
    result = alpha .* f[idx] + beta .* f[idx+1] + gamma .* fpnew[idx] + delta .* fpnew[idx+1]
    return result
end

function interpolate(f, x, order=1)
    # Interpolate positions x (rank order) from field f (rank order)
    # Assume that x are in units of cell widths (dx = 1)
    if order == 0
        interp_fn = interp_zero
    elseif order == 1
        interp_fn = interp_linear
    elseif order == 2
        interp_fn = interp_quad
    elseif order == 3
        interp_fn = interp_cubic
    end
    idx = Int32(div.(x,1))[:]
    dx = mod.(x,1)[:]
    nx = size(x)[1]
    # Left and right derivatives
    fp = grid_deriv(f)
    if ndims(f) == 1
        idiff = (0,1)  # index shift from "idx" to difference
        fi = []
        result = interp_fn(f,fp,idx,idiff,dx)
    elseif ndims(f) == 2
        result = interp_fn(f,fp,idiff,dx)
    elseif ndims(f) == 3
        f000 = f[idx[1],   idx[2]  , idx[3]]
        f010 = f[idx[1],   idx[2]+1, idx[3]]
        f100 = f[idx[1]+1, idx[2]  , idx[3]]
        f110 = f[idx[1]+1, idx[2]+1, idx[3]]
        f001 = f[idx[1],   idx[2]  , idx[3]+1]
        f011 = f[idx[1],   idx[2]+1, idx[3]+1]
        f101 = f[idx[1]+1, idx[2]  , idx[3]+1]
        f111 = f[idx[1]+1, idx[2]+1, idx[3]+1]
        fx00 = interp_fn(f000,f100,dx[1])
        fx01 = interp_fn(f001,f101,dx[1])
        fx10 = interp_fn(f010,f110,dx[1])
        fx11 = interp_fn(f011,f111,dx[1])
        fy0  = interp_fn(fx00,fx10,dx[2])
        fy1  = interp_fn(fx01,fx11,dx[2])
        result = interp_fn(fy0,fy1,dx[3])
    end
    return result
end

#w, mu, th = ang_weights_mu(12, true)