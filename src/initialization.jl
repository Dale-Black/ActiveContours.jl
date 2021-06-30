function ellipsoid_level_set(image_shape::Tuple, center=nothing, semi_axis=nothing)
    if center === nothing
        center = [i ÷ 2 for i in image_shape]
    end

    if semi_axis === nothing
        semi_axis = [i ÷ 2 for i in image_shape]
    end

    if length(center) != length(image_shape)
        error("`center` and `image_shape` must have the same length.")
    end

    if length(semi_axis) != length(image_shape)
        error("`semi_axis` and `image_shape` must have the same length.")
    end

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
        xc, yc, zc = center[1], center[2], center[2]
        rx, ry, rz = semi_axis[1], semi_axis[2], center[3]
        phi = zeros(image_shape)
        for x in 1:image_shape[1]
            for y in 1:image_shape[2]
                for z in 1:image_shape[3]
                    phi[x, y, z] =
                        1 - (((x - xc) / rx)^2 + ((y - yc) / ry)^2 + ((z - zc) / rz)^2)
                end
            end
        end
    else
        error("`image_shape` must be a 2- or 3-tuple.")
    end
    res = phi .> 0
    return res
end
