### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ fac5be42-b57c-4316-b7bb-48e1c65ece44
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate("..")
    using Revise, PlutoUI, PlutoTest, Images, Statistics, ActiveContours
end

# ╔═╡ 10cf431a-d4ac-433f-b5a7-7fb0bcad4bb3
TableOfContents()

# ╔═╡ e242f1df-9865-47fb-95b6-273e51301ada
md"""
# Chan Vese
Julia implementation of [scikit-image's morphological chan-vese segmentation](https://github.com/scikit-image/scikit-image/blob/00177e14097237ef20ed3141ed454bc81b308f82/skimage/segmentation/morphsnakes.py#L214-L314)
"""

# ╔═╡ fd72c0d1-a8b2-42c7-9bbc-b3d6e2dce7d0
md"""
## Helper Functions
"""

# ╔═╡ Cell order:
# ╠═fac5be42-b57c-4316-b7bb-48e1c65ece44
# ╠═10cf431a-d4ac-433f-b5a7-7fb0bcad4bb3
# ╟─e242f1df-9865-47fb-95b6-273e51301ada
# ╟─fd72c0d1-a8b2-42c7-9bbc-b3d6e2dce7d0
