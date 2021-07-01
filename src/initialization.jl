"""
    ellipsoid_level_set(image_shape::Tuple, center=nothing, semi_axis=nothing)

Initialize a 2D or 3D ellipsoid for level set problems. Returns a 
binary array.

# Arguments
- `image_shape::Tuple`: The 2D or 3D shape of the ellipsoid
- `center::Tuple`: Coordinates for the center of the ellipsoid. Defaults to center of 
    of the image.
- `semi-axis::Tuple`: Lengths of the semi-axis of the ellipsoid. Defaults to half of the
    image dimensions.
"""
function ellipsoid_level_set(image_shape::Tuple, center=nothing, semi_axis=nothing)
    if isnothing(center)
        center = [i / 2 for i in image_shape]
    else
        center = collect(center)
    end

    if isnothing(semi_axis)
        semi_axis = [i / 2 for i in image_shape]
    else
        semi_axis = collect(semi_axis)
    end

    length(center) == length(image_shape) ||
        throw(ArgumentError("`center` and `image_shape` must have the same length."))

    length(semi_axis) == length(image_shape) ||
        throw(ArgumentError("`semi_axis` and `image_shape` must have the same length."))

    if length(image_shape) == 2
        xc, yc = center[1], center[2]
        rx, ry = semi_axis[1], semi_axis[2]
        phi = zeros(image_shape)
        for x in 1:image_shape[1]
            for y in 1:image_shape[2]
                phi[x, y] = 1 - ((((x - 1) - xc) / rx)^2 + (((y - 1) - yc) / ry)^2)
            end
        end
    elseif length(image_shape) == 3
        xc, yc, zc = center[1], center[2], center[3]
        rx, ry, rz = semi_axis[1], semi_axis[2], center[3]
        phi = zeros(image_shape)
        for x in 1:image_shape[1]
            for y in 1:image_shape[2]
                for z in 1:image_shape[3]
                    phi[x, y, z] =
                        1 - (
                            (((x - 1) - xc) / rx)^2 +
                            (((y - 1) - yc) / ry)^2 +
                            (((z - 1) - zc) / rz)^2
                        )
                end
            end
        end
    else
        error("`image_shape` must be a 2- or 3-tuple.")
    end
    res = phi .> 0
    return res
end
"""
    lazy_ellipsoid_level_set(image_shape, center=nothing, semi_axis=nothing)

Initialize a 2D or 3D ellipsoid for level set problems. Returns a 
LazySets `Ellipsoid` type. ~~Probably faster than `ellipsoid_level_set`~~.

# Arguments
- `image_shape::Tuple`: The 2D or 3D shape of the ellipsoid
- `center::Tuple`: Coordinates for the center of the ellipsoid. Defaults to center of 
    of the image.
- `semi-axis::Tuple`: Lengths of the semi-axis of the ellipsoid. Defaults to half of the
    image dimensions.
"""
function lazy_ellipsoid_level_set(image_shape::Tuple, center=nothing, semi_axis=nothing)
    if isnothing(center)
        center = [i / 2 for i in image_shape]
    else
        center = collect(center)
    end

    if isnothing(semi_axis)
        semi_axis = [i / 2 for i in image_shape]
    else
        semi_axis = collect(semi_axis)
    end

    length(center) == length(image_shape) ||
        throw(ArgumentError("`center` and `image_shape` must have the same length."))

    length(semi_axis) == length(image_shape) ||
        throw(ArgumentError("`semi_axis` and `image_shape` must have the same length."))

    if length(image_shape) == 2
        rx, ry = semi_axis
        ellipse = Ellipsoid(center, Diagonal([rx^2, ry^2]))
    elseif length(image_shape) == 3
        rx, ry, rz = semi_axis
        ellipse = Ellipsoid(center, Diagonal([rx^2, ry^2, rz^2]))
    else
        error("`image_shape` must be a 2- or 3-tuple.")
    end
    return ellipse
end

"""
    to_array(image_shape::Tuple, E::Ellipsoid)

Turn an `Ellipsoid` type from LazySets.jl into an array
where the pixels within the `Ellipsoid` are 1 and 0 everywhere
else

# Arguments
- `image_shape::Tuple`: The 2D or 3D shape of the `Ellipsoid`
- `E::Ellipsoid`: LazySets `Ellipsoid` type
"""
function to_array(image_shape, E::Ellipsoid)
    B = box_approximation(E)
    Bs = split(B, collect(image_shape))
    res = [center(Bi) in E for Bi in Bs]
    M = reshape(res, image_shape)
    return M
end
