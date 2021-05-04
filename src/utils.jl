# TODO: import this directly from DistanceTransforms.jl
function compute_dtm(img)
    f = ImageMorphology.feature_transform(.!(Bool.(img)))
    foreground_dtm = ImageMorphology.distance_transform(f)
end

function get_curvature(ϕ, idx)
    [dim_y, dim_x] = size(ϕ)
end