"""
1D Bezier interpolation from Hayek et al. (2010) and Kunasz & Auer (1988)
See bezier-h20.jl for more general Bezier interpolation routines that we will mainly use

Calculate control point in a Bezier curve
"""
function bezier_control(s0, sm, sp, dxm, dxp)
    dx_sum = dxm + dxp
    s0prime = (dxp / dx_sum) * (s0 - sm) / dxm + (dxm / dx_sum) * (sp - s0) / dxp
    sc = s0 - 0.5 * dxm * s0prime
    return sc
end

"""
Kunasz & Auer (1988), Equations (7a) - (9c)
Hayek et al. (2010), Appendix A, psim = Psi_d to conform with Kunasz88 notation
"""
    function bezier_interp_coef(s0, sm, sp, dtaum, dtaup, linear=false)
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
        elseif sc < min(sm,s0) || sc > max(sm,s0)
            # if Sc is outside the data range (i.e. overshooting), choose the upwind point as the control point
            psi0 = u2(dtaum) / dtaum^2
            psim = u0(dtaum) - psi0
            psip = 0.0
        end
    end
    return psi0, psim, psip
end

"""
Bezier interpolation for optical depth (Equation A.3 in Hayek+ 2010)
"""
function bezier_interp_tau(chi0, chim, chip, dsm, dsp)
    chic = bezier_control(chi0, chim, chip, dsm, dsp)
    # Limit overshooting
    if max(chim, chi0, chip) == chi0
        chic = chi0
    elseif chic < min(chim, chi0) || chic > max(chim, chi0)
        chic = chim
    end
    ds = dsm + dsp
    # Interpolation
    dtau = ds * (chim + chi0 + chic) / 3
    return dtau
end
