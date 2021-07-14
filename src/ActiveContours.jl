module ActiveContours

using LazySets
using LinearAlgebra

include("./initialization.jl")
include("./morphological_snakes.jl")

export
    # Export initialization.jl functions
    ellipsoid_level_set,
    lazy_ellipsoid_level_set,
    to_array,

    # Export utils.jl functions
    central_diff
end
