### A Pluto.jl notebook ###
# v0.19.18

using Markdown
using InteractiveUtils

# ╔═╡ 0cf5da71-cddb-4143-8bed-184eb98d538f
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("..")
	using Revise, PlutoUI, PlutoTest, Images, Statistics, ActiveContours
end

# ╔═╡ 79db70bb-b291-459a-9b91-9dfabbaee514
TableOfContents()

# ╔═╡ 6ca97834-8b6d-4291-b7b4-1b1c7615328e
md"""
# Chan Vese 
"""

# ╔═╡ e9e8ae1d-bbc9-46a4-a373-ca47963b310f
md"""
## Helper Functions
"""

# ╔═╡ a8d8db5f-2423-4fa4-bdce-cee0171976bf
"""
```julia
_cv_curvature(phi)
```
"""
function _cv_curvature(phi)
    P = padarray(phi, Pad(:replicate, 1, 1))
	fy = (P[2:end, 1:end-1] - P[0:end-2, 1:end-1]) / 2.0
    fx = (P[1:end-1, 2:end] - P[1:end-1, 0:end-2]) / 2.0
    fyy = P[2:end, 1:end-1] + P[0:end-2, 1:end-1] - 2 * phi
    fxx = P[1:end-1, 2:end] + P[1:end-1, 0:end-2] - 2 * phi
    fxy = 0.25 * (P[2:end, 2:end] + P[0:end-2, 0:end-2] - P[0:end-2, 2:end] - P[2:end, 0:end-2])
    grad2 = fx.^2 + fy.^2
	K_num = (fxx .* fy.^2 - 2 * fxy .* fx .* fy + fyy .* fx.^2)
    K_den = (grad2 .* sqrt.(grad2) .+ 1e-8)
	K = K_num ./ K_den
    return K
end

# ╔═╡ 630221d0-9d77-43d8-a772-d25f12cc617e
"""
```julia
_cv_calculate_averages(image, Hphi)
```

Returns the average values 'inside' and 'outside'.
"""
function _cv_calculate_averages(image, Hphi)
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

# ╔═╡ cf823db2-7799-4f48-9179-471b35b78937
"""
```julia
_cv_delta(x, eps=1.0)
```

Returns the result of a regularised dirac delta function of the input value(s).
"""
function _cv_delta(x, eps=1.0)
    eps ./ (eps.^2 .+ x.^2)
end

# ╔═╡ 9d392007-421e-4bea-80ed-c321bac2043c
"""
```julia
_cv_calculate_variation(image, phi, mu, lambda1, lambda2, dt)
```

Returns the variation of level set 'phi' based on algorithm parameters.
"""
function _cv_calculate_variation(image, phi, mu, lambda1, lambda2, dt)
    eta = 1e-16
    P = padarray(phi, Pad(:replicate, 1, 1))
	
    phixp = P[1:end-1, 2:end] - P[1:end-1, 1:end-1]
    phixn = P[1:end-1, 1:end-1] - P[1:end-1, 0:end-2]
    phix0 = (P[1:end-1, 2:end] - P[1:end-1, 0:end-2]) / 2.0

    phiyp = P[2:end, 1:end-1] - P[1:end-1, 1:end-1]
    phiyn = P[1:end-1, 1:end-1] - P[0:end-2, 1:end-1]
    phiy0 = (P[2:end, 1:end-1] - P[0:end-2, 1:end-1]) / 2.0

    C1 = 1 ./ sqrt.(eta .+ phixp.^2 .+ phiy0.^2)
    C2 = 1 ./ sqrt.(eta .+ phixn.^2 .+ phiy0.^2)
    C3 = 1 ./ sqrt.(eta .+ phix0.^2 .+ phiyp.^2)
    C4 = 1 ./ sqrt.(eta .+ phix0.^2 .+ phiyn.^2)

    K = P[1:end-1, 2:end] .* C1 .+ P[1:end-1, 0:end-2] .* C2 .+ P[2:end, 1:end-1] .* C3 .+ P[0:end-2, 1:end-1] .* C4

    Hphi = 1 .* (phi .> 0)
    c1, c2 = _cv_calculate_averages(image, Hphi)

    difference_from_average_term = -1 .* lambda1 .* (image .- c1).^2 .+ lambda2 .* (image .- c2).^2
	
    new_phi = phi + (dt*_cv_delta(phi)) .* (mu*K + difference_from_average_term)
	
	return new_phi ./ (1 .+ mu .* dt .* _cv_delta(phi) .* (C1+C2+C3+C4))
end

# ╔═╡ 7aadad4a-5025-407b-a60f-769b2fda6199
"""
```julia
_cv_heavyside(x, eps=1)
```

Returns the result of a regularised heavyside function of the input value(s).
"""
function _cv_heavyside(x, eps=1)
	return 0.5 .* (1. .+ (2 ./ pi) .* atan.(x ./ eps))
end

# ╔═╡ 4c1907d9-cf8c-448b-945f-9b1afd9572f5
"""
```julia
_cv_difference_from_average_term(image, Hphi, lambda_pos, lambda_neg)
```
"""
function _cv_difference_from_average_term(image, Hphi, lambda_pos, lambda_neg)
    c1, c2 = _cv_calculate_averages(image, Hphi)
    Hinv = 1 .- Hphi
    return (lambda_pos * (image .- c1) .^ 2 .* Hphi + lambda_neg * (image .- c2) .^ 2 .* Hinv)
end

# ╔═╡ 976c59be-f8f0-4e90-8ce8-6c103e8b2416
"""
```julia
_cv_edge_length_term(phi, mu)
```
"""
function _cv_edge_length_term(phi, mu)
    toret = _cv_curvature(phi)
    return mu * toret
end

# ╔═╡ 4f77109c-65e6-4066-8573-2cfd9ecadcc6
"""
```julia
_cv_energy(image, phi, mu, lambda1, lambda2)
```
"""
function _cv_energy(image, phi, mu, lambda1, lambda2)
    H = _cv_heavyside(phi)
    avgenergy = _cv_difference_from_average_term(image, H, lambda1, lambda2)
    lenenergy = _cv_edge_length_term(phi, mu)
    return sum(avgenergy) + sum(lenenergy)
end

# ╔═╡ 4e7b86f3-e77a-4f9f-ae92-e8aec0e61668
md"""
## Chan Vese Function
"""

# ╔═╡ 6b3bc657-f7a5-49f6-a4b5-d6a941905d6d
function chan_vese(image; mu=0.25, lambda1=1, lambda2=1, tol=1e-3, max_iter=500, dt=0.5)
    phi = init_checkerboard(size(image), 5)
    image = image .- minimum(image)
    if maximum(image) != 0
        image = image ./ maximum(image)
    end
    i = 0
    old_energy = _cv_energy(image, phi, mu, lambda1, lambda2)
    energies = []
    phivar = tol + 1
    segmentation = phi .> 0

    while(phivar > tol && i < max_iter)
        # Save old level set values
        oldphi = phi

        # Calculate new level set
        phi = _cv_calculate_variation(image, phi, mu, lambda1, lambda2, dt)
        # phi = _cv_reset_level_set(phi)
        phivar = sqrt(mean((phi-oldphi).^2))

        # Extract energy and compare to previous level set and
        # segmentation to see if continuing is necessary
        segmentation = phi .> 0
        new_energy = _cv_energy(image, phi, mu, lambda1, lambda2)

        # Save old energy values
        push!(energies, old_energy)
        old_energy = new_energy
        i += 1
    end
    return segmentation, phi, energies
end

# ╔═╡ bf4845c4-e92d-4cb2-89b2-d2c7f8ee287d
export chan_vese

# ╔═╡ 161408eb-b994-4dd0-a306-b1fb0d760928
md"""
# Tests
"""

# ╔═╡ 5e82784b-3727-4987-b033-1360c9a73795
test_image = [
  1.0   0.0  -1.0  -0.0   1.0   0.0  -1.0  -0.0   1.0   0.0
  0.0   0.0  -0.0  -0.0   0.0   0.0  -0.0  -0.0   0.0   0.0
 -1.0  -0.0   1.0   0.0  -1.0  -0.0   1.0   0.0  -1.0  -0.0
 -0.0  -0.0   0.0   0.0  -0.0  -0.0   0.0   0.0  -0.0  -0.0
  1.0   0.0  -1.0  -0.0   1.0   0.0  -1.0  -0.0   1.0   0.0
  0.0   0.0  -0.0  -0.0   0.0   0.0  -0.0  -0.0   0.0   0.0
 -1.0  -0.0   1.0   0.0  -1.0  -0.0   1.0   0.0  -1.0  -0.0
 -0.0  -0.0   0.0   0.0  -0.0  -0.0   0.0   0.0  -0.0  -0.0
  1.0   0.0  -1.0  -0.0   1.0   0.0  -1.0  -0.0   1.0   0.0
  0.0   0.0  -0.0  -0.0   0.0   0.0  -0.0  -0.0   0.0   0.0
]

# ╔═╡ 4c2ad612-1587-4669-af61-27c473532739
test_phi = [
  1.0   0.0  -1.0  -0.0   1.0   0.0  -1.0  -0.0   1.0   0.0
  0.0   0.0  -0.0  -0.0   0.0   0.0  -0.0  -0.0   0.0   0.0
 -1.0  -0.0   1.0   0.0  -1.0  -0.0   1.0   0.0  -1.0  -0.0
 -0.0  -0.0   0.0   0.0  -0.0  -0.0   0.0   0.0  -0.0  -0.0
  1.0   0.0  -1.0  -0.0   1.0   0.0  -1.0  -0.0   1.0   0.0
  0.0   0.0  -0.0  -0.0   0.0   0.0  -0.0  -0.0   0.0   0.0
 -1.0  -0.0   1.0   0.0  -1.0  -0.0   1.0   0.0  -1.0  -0.0
 -0.0  -0.0   0.0   0.0  -0.0  -0.0   0.0   0.0  -0.0  -0.0
  1.0   0.0  -1.0  -0.0   1.0   0.0  -1.0  -0.0   1.0   0.0
  0.0   0.0  -0.0  -0.0   0.0   0.0  -0.0  -0.0   0.0   0.0
]

# ╔═╡ f21c1301-0a52-4454-a588-e8bb1195072b
test_Hphi = 1 * (test_phi .> 0)

# ╔═╡ c619f055-f64c-4766-9ca3-11136313121a
begin
	test_mu = 0.25
	test_lambda1 = 1
	test_lambda2 = 1
	test_dt = 0.5
end

# ╔═╡ aed2424b-4e9c-4b1f-b802-9b858558516a
let
	segmentation, phi, energies = chan_vese(test_image; max_iter=200);
	ans = [
		 1  1  1  1  0  0  0  0  0  0
		 1  1  1  1  0  0  0  0  0  0
		 1  1  1  1  0  0  0  0  0  0
		 1  1  1  1  0  0  0  0  0  0
		 0  0  0  0  0  0  0  0  0  0
		 0  0  0  0  0  0  0  0  0  0
		 0  0  0  0  0  0  0  0  0  0
		 0  0  0  0  0  0  0  0  0  0
		 0  0  0  0  0  0  0  0  1  1
		 0  0  0  0  0  0  0  0  1  1
	]
	@test segmentation == ans
end

# ╔═╡ Cell order:
# ╠═0cf5da71-cddb-4143-8bed-184eb98d538f
# ╠═79db70bb-b291-459a-9b91-9dfabbaee514
# ╟─6ca97834-8b6d-4291-b7b4-1b1c7615328e
# ╟─e9e8ae1d-bbc9-46a4-a373-ca47963b310f
# ╠═a8d8db5f-2423-4fa4-bdce-cee0171976bf
# ╠═630221d0-9d77-43d8-a772-d25f12cc617e
# ╠═9d392007-421e-4bea-80ed-c321bac2043c
# ╠═cf823db2-7799-4f48-9179-471b35b78937
# ╠═7aadad4a-5025-407b-a60f-769b2fda6199
# ╠═4c1907d9-cf8c-448b-945f-9b1afd9572f5
# ╠═976c59be-f8f0-4e90-8ce8-6c103e8b2416
# ╠═4f77109c-65e6-4066-8573-2cfd9ecadcc6
# ╟─4e7b86f3-e77a-4f9f-ae92-e8aec0e61668
# ╠═6b3bc657-f7a5-49f6-a4b5-d6a941905d6d
# ╠═bf4845c4-e92d-4cb2-89b2-d2c7f8ee287d
# ╟─161408eb-b994-4dd0-a306-b1fb0d760928
# ╠═5e82784b-3727-4987-b033-1360c9a73795
# ╠═4c2ad612-1587-4669-af61-27c473532739
# ╠═f21c1301-0a52-4454-a588-e8bb1195072b
# ╠═c619f055-f64c-4766-9ca3-11136313121a
# ╠═aed2424b-4e9c-4b1f-b802-9b858558516a
