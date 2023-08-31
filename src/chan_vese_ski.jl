"""
    chan_vese_ski(image, phi; mu=0.25, lambda1=1, lambda2=1, tol=1e-3, max_iter=500, dt=0.5)

Julia implementation of [scikit-image's chan-vese segmentation](https://github.com/scikit-image/scikit-image/blob/v0.19.2/skimage/segmentation/_chan_vese.py#L175-L347)
...
"""
function chan_vese_ski(image, phi; mu=0.25, lambda1=1, lambda2=1, tol=1e-3, max_iter=500, dt=0.5)
    image = image .- minimum(image)
    if maximum(image) != 0
        image = image ./ maximum(image)
    end
    i = 0
    old_energy = energy(image, phi, mu, lambda1, lambda2)
    energies = []
    phivar = tol + 1
    segmentation = phi .> 0

    while (phivar > tol && i < max_iter)
        # Save old level set values
        oldphi = phi

        # Calculate new level set
        phi = calculate_variation(image, phi, mu, lambda1, lambda2, dt)
        # phi = _cv_reset_level_set(phi)
        phivar = sqrt(mean((phi - oldphi) .^ 2))

        # Extract energy and compare to previous level set and
        # segmentation to see if continuing is necessary
        segmentation = phi .> 0
        new_energy = energy(image, phi, mu, lambda1, lambda2)

        # Save old energy values
        push!(energies, old_energy)
        old_energy = new_energy
        i += 1
    end
    return segmentation, phi, energies
end

export chan_vese_ski

"""
    curvature(phi)
"""
function curvature(phi)
    P = padarray(phi, Pad(:replicate, 1, 1))
    fy = (P[2:end, 1:end-1] - P[0:end-2, 1:end-1]) / 2.0
    fx = (P[1:end-1, 2:end] - P[1:end-1, 0:end-2]) / 2.0
    fyy = P[2:end, 1:end-1] + P[0:end-2, 1:end-1] - 2 * phi
    fxx = P[1:end-1, 2:end] + P[1:end-1, 0:end-2] - 2 * phi
    fxy = 0.25 * (P[2:end, 2:end] + P[0:end-2, 0:end-2] - P[0:end-2, 2:end] - P[2:end, 0:end-2])
    grad2 = fx .^ 2 + fy .^ 2
    K = (fxx .* fy .^ 2 - 2 * fxy .* fx .* fy + fyy .* fx .^ 2) ./ (grad2 .* sqrt.(grad2) .+ 1e-8)
    return K
end

export curvature

"""
    calculate_averages(image, Hphi)

Returns the average values 'inside' and 'outside'.
"""
function calculate_averages(image, Hphi)
    H = Hphi
    Hinv = 1 .- H
    Hsum = sum(H)
    Hinvsum = sum(Hinv)
    avg_inside = sum(image .* H)
    avg_outside = sum(image .* Hinv)
    if Hsum != 0
        avg_inside /= Hsum
    end
    if Hinvsum != 0
        avg_outside /= Hinvsum
    end
    return (avg_inside, avg_outside)
end

"""
    delta(x; eps=1)

Returns the result of a regularised Dirac delta function of the input value(s).
"""
function delta(x; eps=1)
    eps ./ (eps .^ 2 .+ x .^ 2)
end

"""
    calculate_variation(image, phi, mu, lambda1, lambda2, dt)

Returns the variation of level set 'phi' based on algorithm parameters.
"""
function calculate_variation(image, phi, mu, lambda1, lambda2, dt; eta = 1e-16)
    P = padarray(phi, Pad(:replicate, 1, 1))

    phixp = P[1:end-1, 2:end] - P[1:end-1, 1:end-1]
    phixn = P[1:end-1, 1:end-1] - P[1:end-1, 0:end-2]
    phix0 = (P[1:end-1, 2:end] - P[1:end-1, 0:end-2]) / 2.0

    phiyp = P[2:end, 1:end-1] - P[1:end-1, 1:end-1]
    phiyn = P[1:end-1, 1:end-1] - P[0:end-2, 1:end-1]
    phiy0 = (P[2:end, 1:end-1] - P[0:end-2, 1:end-1]) / 2.0

    C1 = 1 ./ sqrt.(eta .+ phixp .^ 2 .+ phiy0 .^ 2)
    C2 = 1 ./ sqrt.(eta .+ phixn .^ 2 .+ phiy0 .^ 2)
    C3 = 1 ./ sqrt.(eta .+ phix0 .^ 2 .+ phiyp .^ 2)
    C4 = 1 ./ sqrt.(eta .+ phix0 .^ 2 .+ phiyn .^ 2)

    K = P[1:end-1, 2:end] .* C1 .+ P[1:end-1, 0:end-2] .* C2 .+ P[2:end, 1:end-1] .* C3 .+ P[0:end-2, 1:end-1] .* C4

    Hphi = 1 .* (phi .> 0)
    c1, c2 = calculate_averages(image, Hphi)

    difference_from_average_term = -1 .* lambda1 .* (image .- c1) .^ 2 .+ lambda2 .* (image .- c2) .^ 2

    new_phi = phi + (dt * delta(phi)) .* (mu * K + difference_from_average_term)

    return new_phi ./ (1 .+ mu .* dt .* delta(phi) .* (C1 + C2 + C3 + C4))
end

"""
    heavyside(x; eps=1)

Returns the result of a regularised heavyside function of the input value(s).
"""
function heavyside(x; eps=1)
    return 0.5 .* (1.0 .+ (2 ./ pi) .* atan.(x ./ eps))
end

"""
    difference_from_average(image, Hphi, lambda_pos, lambda_neg)
"""
function difference_from_average(image, Hphi, lambda_pos, lambda_neg)
    c1, c2 = calculate_averages(image, Hphi)
    Hinv = 1 .- Hphi
    return (lambda_pos * (image .- c1) .^ 2 .* Hphi + lambda_neg * (image .- c2) .^ 2 .* Hinv)
end

"""
    edge_length(phi, mu)
"""
function edge_length(phi, mu)
    toret = curvature(phi)
    return mu * toret
end

"""
    energy(image, phi, mu, lambda1, lambda2)
"""
function energy(image, phi, mu, lambda1, lambda2)
    H = heavyside(phi)
    avgenergy = difference_from_average(image, H, lambda1, lambda2)
    lenenergy = edge_length(phi, mu)
    return sum(avgenergy) + sum(lenenergy)
end

