### A Pluto.jl notebook ###
# v0.19.18

using Markdown
using InteractiveUtils

# ╔═╡ a6b0911a-9db4-11ed-2e08-e10ba1344306
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("..")
	using Revise, PlutoUI, PlutoTest
end

# ╔═╡ 112cb9ba-6dd0-4881-a67e-b104dc1a38c5
TableOfContents()

# ╔═╡ 18da3cf3-c3a5-4d43-a584-dc0b63cf7a8e
function init_checkerboard(image_size, square_size)
	yv = reshape(collect(range(1, image_size[1])), (image_size[1], 1))
    xv = collect(range(1, image_size[2]))
    sf = pi / square_size
    xv = xv .* sf
    yv = yv .* sf
    return sin.(yv) .* sin.(xv)'
end

# ╔═╡ 00ee61b7-bd99-4a01-a985-af6772f02a95
export init_checkerboard

# ╔═╡ 7c608bbf-0105-4036-8e37-6adddc0fc7d6
round.(init_checkerboard((10, 10), 2))

# ╔═╡ Cell order:
# ╠═a6b0911a-9db4-11ed-2e08-e10ba1344306
# ╠═112cb9ba-6dd0-4881-a67e-b104dc1a38c5
# ╠═18da3cf3-c3a5-4d43-a584-dc0b63cf7a8e
# ╠═00ee61b7-bd99-4a01-a985-af6772f02a95
# ╠═7c608bbf-0105-4036-8e37-6adddc0fc7d6
