### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 287234e2-3231-11ee-326c-211056a128fc
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(pwd())
	using Revise
	using ActiveContours, TestImages, CairoMakie, PlutoUI, DICOM, Images, ImageSegmentation, Statistics
	using StatsBase: countmap
end

# ╔═╡ 888681c9-afa8-45ba-bc04-fcf9b0180e1e
TableOfContents()

# ╔═╡ d1cc343d-0278-4c2a-baed-26387bbbb645
md"""
# Cameraman
"""

# ╔═╡ 3fb0b020-3bcb-4090-818e-82feaf187f14
cam_img = testimage("cameraman")

# ╔═╡ 8e4ba501-2fe4-4b2c-8d43-2d9eb5422845
cam_mask = chan_vese(cam_img);

# ╔═╡ 701b69ea-ca00-42d5-86c9-4cd8420f9db0
heatmap(rotr90(cam_mask), colormap = :grays)

# ╔═╡ 817cfbb5-de41-4ddc-a008-f2f680418b9b
# chan_vese_ski(cam_img, init_checkerboard(size(cam_img), 5); max_iter=200)

# ╔═╡ 3ffbba76-f898-42c2-9559-05c9505b264c
md"""
# Brain
"""

# ╔═╡ a0faee44-4089-4a6a-a0b5-8f1e3fa0a7a8
brain_img = testimage("mri-stack");

# ╔═╡ d5c792aa-1f9b-42e8-8432-0d1d487ad377
brain_mask = chan_vese(brain_img; max_iter = 500);

# ╔═╡ 2d80ae14-1b50-47a1-8315-3e13bfdd8080
@bind z1 PlutoUI.Slider(axes(brain_img, 3), show_value = true)

# ╔═╡ b594d922-2c6f-41f5-95af-55693cf7b85d
let
	f = Figure(resolution = (1200, 800))
	ax = CairoMakie.Axis(
		f[1, 1],
		title = "Raw Brain MRI"
	)
	heatmap!(rotr90(brain_img[:, :, z1]))

	ax = CairoMakie.Axis(
		f[1, 2],
		title = "Mask Brain MRI"
	)
	heatmap!(rotr90(brain_mask[:, :, z1]), colormap = :grays)
	
	f
end

# ╔═╡ 17de12c2-b54f-445f-9742-0cce060f04a5
md"""
# DICOM
"""

# ╔═╡ b7ab28bf-115a-4055-ba0b-2091a7db66a0
dir = "/Users/daleblack/Desktop/dcms/2"

# ╔═╡ b8d8df5c-1203-4f1b-9a98-4992fab4b6c1
dcm = dcmdir_parse(dir);

# ╔═╡ ff536518-23a6-4016-8342-accea257e3be
function load_dcm_array(dcm_data::Vector{DICOM.DICOMData})
    return array = cat(
        [dcm_data[i][(0x7fe0, 0x0010)] for i in 1:length(dcm_data)]...; dims=3
    )
end

# ╔═╡ 09af52e5-fc6f-4c06-ac1d-0cf0ac790d4f
heart_img = load_dcm_array(dcm);

# ╔═╡ e3c5118a-9121-4c92-b28c-29c1e4ec30f9
heart_mask = chan_vese(heart_img);

# ╔═╡ 2571f84c-6f79-4a23-ade2-0b0d8d062e63
@bind z2 PlutoUI.Slider(axes(heart_mask, 3), show_value = true)

# ╔═╡ 4f1db59c-8d89-453b-8981-281a66201f8b
let
	f = Figure(resolution = (2000, 1400))
	ax = CairoMakie.Axis(
		f[1, 1],
		title = "Raw Heart CT"
	)
	heatmap!(rotr90(heart_img[:, :, z2]), colormap = :grays)

	ax = CairoMakie.Axis(
		f[1, 2],
		title = "Mask Heart CT"
	)
	heatmap!(rotr90(heart_mask[:, :, z2]), colormap = :grays)
	
	f
end

# ╔═╡ 9779e00d-238f-4688-b2f0-d99a2633958a
function centroids_from_mask(mask)
	cc_labels = label_components(mask)
	largest_connected_component, _ = sort(collect(pairs(countmap(cc_labels[cc_labels .!= 0]))), by=x->x[2], rev=true)
	largest_connected_indices = findall(cc_labels .== largest_connected_component[1])

	new_mask = zeros(size(mask))
	for i in largest_connected_indices
		new_mask[i] = 1
	end
	centroids = Int.(round.(component_centroids(label_components(new_mask))[end]))
end

# ╔═╡ 84375483-0aa0-4059-8dea-95fbce3cb1c3
centroids = centroids_from_mask(heart_mask)

# ╔═╡ ca8e7273-59c5-4a64-8e78-fada39d397b4
let
	f = Figure()
	ax = CairoMakie.Axis(
		f[1, 1]
	)
	heatmap!(rotr90(heart_img[:, :, 160]), colormap = :grays)
	scatter!([centroids[2]], [centroids[1]], markersize = 10, color = :red)
	f
end

# ╔═╡ 2418b23a-8e8a-42bf-941b-5c5d0898e268
function initial_level_set(shape::Tuple{Int64, Int64, Int64})
    x₀ = reshape(collect(0:shape[begin]-1), shape[begin], 1, 1)
    y₀ = reshape(collect(0:shape[begin+1]-1), 1, shape[begin+1], 1)
    z₀ = reshape(collect(0:shape[begin+2]-1), 1, 1, shape[begin+2])
    𝚽₀ = @. sin(pi / 5 * x₀) * sin(pi / 5 * y₀) * sin(pi / 5 * z₀)
end

# ╔═╡ 92166e22-ed47-4f93-be18-fc126b40960e
level_set2 = initial_level_set(size(heart_img)) .* heart_mask;

# ╔═╡ 2e33f4a2-af87-401e-b6cf-826d872bc6bc
heatmap(level_set2[:, :, 100])

# ╔═╡ d004d912-ac1b-4579-90d7-04485713479b
heart_mask2 = chan_vese(heart_img; init_level_set = level_set2);

# ╔═╡ 7e7c530c-7ae9-49bd-8f78-f2a0b8ebcac3
heatmap(heart_mask2[:, :, 100])

# ╔═╡ Cell order:
# ╠═287234e2-3231-11ee-326c-211056a128fc
# ╠═888681c9-afa8-45ba-bc04-fcf9b0180e1e
# ╟─d1cc343d-0278-4c2a-baed-26387bbbb645
# ╠═3fb0b020-3bcb-4090-818e-82feaf187f14
# ╠═8e4ba501-2fe4-4b2c-8d43-2d9eb5422845
# ╠═701b69ea-ca00-42d5-86c9-4cd8420f9db0
# ╠═817cfbb5-de41-4ddc-a008-f2f680418b9b
# ╟─3ffbba76-f898-42c2-9559-05c9505b264c
# ╠═a0faee44-4089-4a6a-a0b5-8f1e3fa0a7a8
# ╠═d5c792aa-1f9b-42e8-8432-0d1d487ad377
# ╟─2d80ae14-1b50-47a1-8315-3e13bfdd8080
# ╟─b594d922-2c6f-41f5-95af-55693cf7b85d
# ╟─17de12c2-b54f-445f-9742-0cce060f04a5
# ╠═b7ab28bf-115a-4055-ba0b-2091a7db66a0
# ╠═b8d8df5c-1203-4f1b-9a98-4992fab4b6c1
# ╠═ff536518-23a6-4016-8342-accea257e3be
# ╠═09af52e5-fc6f-4c06-ac1d-0cf0ac790d4f
# ╠═e3c5118a-9121-4c92-b28c-29c1e4ec30f9
# ╟─2571f84c-6f79-4a23-ade2-0b0d8d062e63
# ╟─4f1db59c-8d89-453b-8981-281a66201f8b
# ╠═9779e00d-238f-4688-b2f0-d99a2633958a
# ╠═84375483-0aa0-4059-8dea-95fbce3cb1c3
# ╠═ca8e7273-59c5-4a64-8e78-fada39d397b4
# ╠═92166e22-ed47-4f93-be18-fc126b40960e
# ╠═2e33f4a2-af87-401e-b6cf-826d872bc6bc
# ╠═2418b23a-8e8a-42bf-941b-5c5d0898e268
# ╠═d004d912-ac1b-4579-90d7-04485713479b
# ╠═7e7c530c-7ae9-49bd-8f78-f2a0b8ebcac3
