module ActiveContours

using LazySets

include("./initialization.jl")
include("./morphological_snakes.jl")

export
    # Export initialization.jl functions
    ellipsoid_level_set,
    lazy_ellipsoid_level_set,
    to_array
end
