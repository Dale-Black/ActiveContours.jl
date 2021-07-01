include("./imports.jl")

@testset ExtendedTestSet "ellipsoid_level_set" begin
    @testset ExtendedTestSet "ellipsoid_level_set" begin
        image_shape = (10, 10)
        answer2D = [
            0 0 0 0 0 0 0 0 0 0
            0 0 0 1 1 1 1 1 0 0
            0 0 1 1 1 1 1 1 1 0
            0 1 1 1 1 1 1 1 1 1
            0 1 1 1 1 1 1 1 1 1
            0 1 1 1 1 1 1 1 1 1
            0 1 1 1 1 1 1 1 1 1
            0 1 1 1 1 1 1 1 1 1
            0 0 1 1 1 1 1 1 1 0
            0 0 0 1 1 1 1 1 0 0
        ]
        test = ellipsoid_level_set(image_shape)
        @test test == answer2D
    end

    @testset ExtendedTestSet "ellipsoid_level_set" begin
        image_shape = (10, 10, 2)
        a1 = [
            0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0
        ]

        a2 = [
            0 0 0 0 0 0 0 0 0 0
            0 0 0 1 1 1 1 1 0 0
            0 0 1 1 1 1 1 1 1 0
            0 1 1 1 1 1 1 1 1 1
            0 1 1 1 1 1 1 1 1 1
            0 1 1 1 1 1 1 1 1 1
            0 1 1 1 1 1 1 1 1 1
            0 1 1 1 1 1 1 1 1 1
            0 0 1 1 1 1 1 1 1 0
            0 0 0 1 1 1 1 1 0 0
        ]
        answer3D = cat(a1, a2; dims=3)
        test = ellipsoid_level_set(image_shape)
        @test test == answer3D
    end
end
