### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 68b85236-14eb-42ff-bbb7-3e21807ecec0
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using Revise, Images, TestImages, PlutoUI, CairoMakie, ActiveContours
end

# ╔═╡ 6a1bb0eb-6a7c-442a-a386-7f7ea1640f0b
TableOfContents()

# ╔═╡ 598c7ac2-b30e-46c0-ae48-9a88976fd691
md"""
## Load image
"""

# ╔═╡ e7f68608-2fff-418c-87fa-a4d1272fef07
begin
	# img = load(joinpath(dirname(pwd()), "images/blobs.png"))
	img = TestImages.shepp_logan(256)
	img = convert(Array{Float64}, Gray.(img))
end

# ╔═╡ f734374a-7227-43b6-8aa5-ed1372a4c45e
phi0 = init_checkerboard(size(img), 5)

# ╔═╡ 117afcfd-6c59-4b07-a5df-b56080c3a31a
md"""
## Add functions
"""

# ╔═╡ ec54e56d-a543-413b-8fbb-ce6eedf22edb
segmentation, phi, energies = chan_vese(img, phi0; max_iter=200);

# ╔═╡ 69d20de2-da38-4dfe-bc2e-d765cf27f66b
md"""
## Visualize
"""

# ╔═╡ cc683ab4-0d3b-4cbc-914b-0545db0a84a8
let
	f = Figure()

	ax1 = CairoMakie.Axis(
		f[1, 1],
		title="Original Image"
	)
	heatmap!(img, colormap=:grays)

	ax2 = CairoMakie.Axis(
		f[1, 2],
		title="Final Level Set"
	)
	heatmap!(phi, colormap=:grays)

	ax3 = CairoMakie.Axis(
		f[2, 1],
		title="Chan-Vese Segmentation"
	)
	heatmap!(segmentation, colormap=:grays)

	ax4 = CairoMakie.Axis(
		f[2, 2],
		title="Evolution of energy"
	)
	lines!(Float32.(energies))

	f
end

# ╔═╡ Cell order:
# ╠═68b85236-14eb-42ff-bbb7-3e21807ecec0
# ╠═6a1bb0eb-6a7c-442a-a386-7f7ea1640f0b
# ╟─598c7ac2-b30e-46c0-ae48-9a88976fd691
# ╠═e7f68608-2fff-418c-87fa-a4d1272fef07
# ╠═f734374a-7227-43b6-8aa5-ed1372a4c45e
# ╟─117afcfd-6c59-4b07-a5df-b56080c3a31a
# ╠═ec54e56d-a543-413b-8fbb-ce6eedf22edb
# ╟─69d20de2-da38-4dfe-bc2e-d765cf27f66b
# ╟─cc683ab4-0d3b-4cbc-914b-0545db0a84a8
