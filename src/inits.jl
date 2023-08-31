"""
	init_checkerboard(image_size, square_size)

Returns a boolean checkerboard image of size `image_size` with squares of size `square_size`.

"""
function init_checkerboard(image_size, square_size)
	rows, cols = image_size
	sf = pi / square_size
	
	yv = sin.(collect(1:rows)' .* sf)
	xv = sin.(collect(1:cols) .* sf)
    return yv .* xv
end

export init_checkerboard