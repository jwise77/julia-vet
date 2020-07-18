# Calculate Gaussian angular quadrature: method of discrete ordinates
# (Bruls et al. 1999, Appendix B)
#
# mu = cos(theta)
#
# nmu=2: point for each quadrant.
# Athena-RT (Davis et al. 2012) uses up to nmu=12
# Results in a na = nmu*(nmu+2) total number of rays

using Printf

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
    # Hayek et al. (2010), Appendix A
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
        else if sc < min(sm,s0) && sc > max(sm,s0)
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
    else if chic < min(chim, chi0) && chic > max(chim, chi0)
        chic = chim
    end
    ds = dsm + dsp
    # Interpolation
    dtau = ds * (chim + chi0 + chic) / 3
    return dtau
end

#w, mu, th = ang_weights_mu(12, true)