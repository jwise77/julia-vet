using Printf

"""
Calculate Gaussian angular quadrature: method of discrete ordinates
(Bruls et al. 1999, Appendix B)

Returns (weights, mu, theta)

mu = cos(theta)

nmu=2: point for each quadrant.
Athena-RT (Davis et al. 2012) uses up to nmu=12
Results in a na = nmu*(nmu+2) total number of rays
"""
function ang_weights_mu(nmu, north_only=true, verbose=false)

    if nmu == 2
        # Trivial answer for one point
        mu = sqrt(1/3)
        w = 1.0
        theta = acos(mu[1])
        return [w,w], [mu,-mu], [theta,-theta]
    elseif nmu % 2 == 1
        println("Number of polar angles must be even.  Setting n_mu = $nmu -> $(nmu+1)")
        nmu += 1
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
    if north_only
        w = reverse(wN)
        mu = reverse(muN)
        theta = reverse(thetaN)
    else
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
    end

    if verbose
        println("mu = $mu")
        println("theta = $theta")
        println("w = $w")
        @printf "sum(w) = %f\n" sum(w)
    end

    return w, mu, theta
    
end

"""
Calculate the angle coordinates (mu_i, mu_j, mu_k) and quadrature weights

Input: nmu (number of polar angles, theta)
Output: ray_info dictionary

See Bruls+ (1999), Appendix B. 
Note (Aug 2021): For testing, we take the simplification of having z as the preferred direction because we're on a Cartesian grid. This allows us to neglect the rotation invariance and define the points on co-latitude (at theta) circle evenly in longitude (phi).  In the future, we need to solve a system of linear equations for the class weights.
"""
function calculate_ray_info(nmu::Int)
    if nmu <= 1
        error("nmu must be greater than 1")
    elseif nmu % 2 == 1
        println("WARNING: Number of polar angles must be even.  Setting n_mu = $nmu -> $(nmu+1)")
        nmu += 1
    end

    na = nmu * (nmu+2)
    w_z, mu_z, theta = ang_weights_mu(nmu)
    
    # Ray angle and weight arrays
    mu = zeros(na, 3)
    w = zeros(na)

    # Assign phi. Work on north hemisphere and copy to southern
    # n_j (= 2*nmu + 4 - 4*j) points on the j-th circle in the unit sphere
    istart = 0
    iend = 0
    for j = 1:div(nmu,2)
        nj = Int(2*nmu + 4 - 4*j)
        # Update array bounds for this co-latitude circle
        istart = iend+1
        iend += nj
        dphi = 2*pi/nj
        phi = (collect(1:nj) .- 0.5) .* dphi
        mu[istart:iend,1] = sin(theta[j]) .* cos.(phi)
        mu[istart:iend,2] = sin(theta[j]) .* sin.(phi)
        mu[istart:iend,3] = ones(nj) .* mu_z[j]
        w[istart:iend] = ones(nj) .* w_z[j] ./ nj
    end
    # Copy to southern hemisphere
    mu[iend+1:end,1] = mu[1:iend,1]
    mu[iend+1:end,2] = mu[1:iend,2]
    mu[iend+1:end,3] = -mu[1:iend,3]
    w[iend+1:end] = w[1:iend]

    # Calculate intersection with cell faces (take cell width ds = 1)
    # Distance to neighboring xy-, xz-, yz-planes.
    # The minimum is the ray segment length.
    ds_xy = sign.(mu[:,3]) ./ mu[:,3]
    ds_xz = sign.(mu[:,2]) ./ mu[:,2]
    ds_yz = sign.(mu[:,1]) ./ mu[:,1]
    ds = min.(ds_xy, ds_xz, ds_yz)
    dr = ds .* mu

    # Cartesian direction (i.e. x,y,z) to the nearest plane
    idir = dropdims(map(a -> a[2], argmax(abs.(dr), dims=2)), dims=2)
    # Pointing in positive or negative direction
    isign = [sign(dr[i,idir[i]]) for i = 1:na]

    ray_info = Dict("nmu"=>nmu, "na"=>nmu*(nmu+2), "w"=>w, "mu"=>mu, 
        "ds"=>ds, "dr"=>dr, "idir"=>idir, "isign"=>isign)

    return ray_info
end
