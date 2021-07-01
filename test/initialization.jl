include("./imports.jl")

@testset ExtendedTestSet "ellipsoid_level_set" begin
    @testset ExtendedTestSet "ellipsoid_level_set" begin
        image_shape = (10, 10)
        answer = [
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
        @test test == answer
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
        answer = cat(a1, a2; dims=3)
        test = ellipsoid_level_set(image_shape)
        @test test == answer
    end
end

@testset ExtendedTestSet "lazy_ellipsoid_level_set" begin
    @testset ExtendedTestSet "lazy_ellipsoid_level_set" begin
        image_shape = (10, 10)
        answer = LazySets.Ellipsoid([5.0, 5.0], Diagonal([25.0, 25.0]))
        test = lazy_ellipsoid_level_set(image_shape)
        @test test == answer
    end

    @testset ExtendedTestSet "lazy_ellipsoid_level_set" begin
        image_shape = (10, 10, 2)
        # a1 = [
        #     0 0 0 0 0 0 0 0 0 0
        #     0 0 1 1 1 1 1 1 0 0
        #     0 1 1 1 1 1 1 1 1 0
        #     0 1 1 1 1 1 1 1 1 0
        #     0 1 1 1 1 1 1 1 1 0
        #     0 1 1 1 1 1 1 1 1 0
        #     0 1 1 1 1 1 1 1 1 0
        #     0 1 1 1 1 1 1 1 1 0
        #     0 0 1 1 1 1 1 1 0 0
        #     0 0 0 0 0 0 0 0 0 0
        # ]

        # a2 = [
        #     0 0 0 0 0 0 0 0 0 0
        #     0 0 1 1 1 1 1 1 0 0
        #     0 1 1 1 1 1 1 1 1 0
        #     0 1 1 1 1 1 1 1 1 0
        #     0 1 1 1 1 1 1 1 1 0
        #     0 1 1 1 1 1 1 1 1 0
        #     0 1 1 1 1 1 1 1 1 0
        #     0 1 1 1 1 1 1 1 1 0
        #     0 0 1 1 1 1 1 1 0 0
        #     0 0 0 0 0 0 0 0 0 0
        # ]
        answer = LazySets.Ellipsoid([5.0, 5.0, 1.0], Diagonal([25.0, 25.0, 1.0]))
        test = lazy_ellipsoid_level_set(image_shape)
        @test test == answer
    end
end

@testset ExtendedTestSet "to_array" begin
    @testset ExtendedTestSet "to_array" begin
        image_shape = (10, 10)
        answer = [
            0 0 0 1 1 1 1 0 0 0
            0 1 1 1 1 1 1 1 1 0
            0 1 1 1 1 1 1 1 1 0
            1 1 1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1 1 1
            0 1 1 1 1 1 1 1 1 0
            0 1 1 1 1 1 1 1 1 0
            0 0 0 1 1 1 1 0 0 0
        ]
        ellipse = LazySets.Ellipsoid([5.0, 5.0], Diagonal([25.0, 25.0]))
        test = to_array(image_shape, ellipse)
        @test test == answer
    end

    @testset ExtendedTestSet "lazy_ellipsoid_level_set" begin
        image_shape = (10, 10, 2)
        a1 = [
            0 0 0 0 0 0 0 0 0 0
            0 0 1 1 1 1 1 1 0 0
            0 1 1 1 1 1 1 1 1 0
            0 1 1 1 1 1 1 1 1 0
            0 1 1 1 1 1 1 1 1 0
            0 1 1 1 1 1 1 1 1 0
            0 1 1 1 1 1 1 1 1 0
            0 1 1 1 1 1 1 1 1 0
            0 0 1 1 1 1 1 1 0 0
            0 0 0 0 0 0 0 0 0 0
        ]

        a2 = [
            0 0 0 0 0 0 0 0 0 0
            0 0 1 1 1 1 1 1 0 0
            0 1 1 1 1 1 1 1 1 0
            0 1 1 1 1 1 1 1 1 0
            0 1 1 1 1 1 1 1 1 0
            0 1 1 1 1 1 1 1 1 0
            0 1 1 1 1 1 1 1 1 0
            0 1 1 1 1 1 1 1 1 0
            0 0 1 1 1 1 1 1 0 0
            0 0 0 0 0 0 0 0 0 0
        ]
        answer = cat(a1, a2, dims=3)
        ellipse = LazySets.Ellipsoid([5.0, 5.0, 1.0], Diagonal([25.0, 25.0, 1.0]))
        test = to_array(image_shape, ellipse)
        @test test == answer
    end
end
