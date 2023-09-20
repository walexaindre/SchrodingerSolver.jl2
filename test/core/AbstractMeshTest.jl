
@testset "2D AbstractMesh" begin
    @testset "Basic Functionality I" begin
        for x_size in 7:11
            for y_size in 6:11
                mesh2d = MetaMesh2D(x_size, y_size)

                @test linear_size(mesh2d) == y_size * x_size

                @test step_x(1, mesh2d) == (y_size, 2)
                @test step_x(y_size, mesh2d) == (y_size - 1, 1)
                @test step_x(3, mesh2d) == (2, 4)
                @test step_y(1, mesh2d) == (x_size, 2)
                @test step_y(x_size, mesh2d) == (x_size - 1, 1)
                @test step_y(3, mesh2d) == (2, 4)

                @test linear_indexing(3, 1, mesh2d) == 3
                @test linear_indexing(3, 2, mesh2d) == x_size + 3

                @test linear_indexing(x_size * y_size, mesh2d) == (x_size, y_size)
                @test linear_indexing(x_size, y_size, mesh2d) == x_size * y_size

                @test linear_indexing(x_size + 3, mesh2d) == (3, 2)
                @test linear_indexing(2 * x_size + 5, mesh2d) == (5, 3)

                @test step_x(1, 2, mesh2d) == (y_size - 1, 3)
                @test step_x(3, 2, mesh2d) == (1, 5)
                @test step_x(4, 2, mesh2d) == (2, 6)

                @test step_y(x_size, 3, mesh2d) == (x_size - 3, 3)
                @test step_y(x_size - 3, 3, mesh2d) == (x_size - 6, x_size)
                @test step_y(4, 2, mesh2d) == (2, 6)
            end
        end
    end
end
