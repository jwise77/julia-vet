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
