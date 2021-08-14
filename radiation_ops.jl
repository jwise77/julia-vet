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

