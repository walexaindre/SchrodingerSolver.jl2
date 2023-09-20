@testset "1D AbstractMesh" begin
    @testset "Basic Functionality I" begin
        for x_size in 9:14
            mesh1d = MetaMesh1D(x_size)

            @test linear_size(mesh1d) == x_size

            #Step x direction
            @test step_x(1, mesh1d) == (x_size, 2)
            @test step_x(x_size, mesh1d) == (x_size - 1, 1)
            @test step_x(3, mesh1d) == (2, 4)

            #Linear indexing back and forth
            @test linear_indexing(1, mesh1d) == 1
            @test linear_indexing(6, mesh1d) == 6

            @test linear_indexing(x_size,mesh1d) == x_size
            @test linear_indexing(x_size-2,mesh1d) == x_size-2

            #Step x direction length more than 1
            @test step_x(1, 2, mesh1d) == (x_size - 1, 3)
            @test step_x(3, 2, mesh1d) == (1, 5)
            @test step_x(4, 2, mesh1d) == (2, 6)
            @test step_x(x_size,2,mesh1d) == (x_size-2,2)

            @test step_x(4,3,mesh1d) == (1,7)
            @test step_x(x_size,3,mesh1d) == (x_size-3,3)
        end
    end
end

@testset "2D AbstractMesh" begin
    @testset "Basic Functionality I" begin
        for x_size in 7:11
            for y_size in 6:11
                mesh2d = MetaMesh2D(x_size, y_size)

                @test linear_size(mesh2d) == y_size * x_size

                #Step x direction
                @test step_x(1, mesh2d) == (y_size, 2)
                @test step_x(y_size, mesh2d) == (y_size - 1, 1)
                @test step_x(3, mesh2d) == (2, 4)

                #Step y direction
                @test step_y(1, mesh2d) == (x_size, 2)
                @test step_y(x_size, mesh2d) == (x_size - 1, 1)
                @test step_y(3, mesh2d) == (2, 4)

                #Linear indexing back and forth
                @test linear_indexing(3, 1, mesh2d) == 3
                @test linear_indexing(3, 2, mesh2d) == x_size + 3

                @test linear_indexing(x_size * y_size, mesh2d) == (x_size, y_size)
                @test linear_indexing(x_size, y_size, mesh2d) == x_size * y_size

                @test linear_indexing(x_size + 3, mesh2d) == (3, 2)
                @test linear_indexing(2 * x_size + 5, mesh2d) == (5, 3)

                #Step x direction length more than 1
                @test step_x(1, 2, mesh2d) == (y_size - 1, 3)
                @test step_x(3, 2, mesh2d) == (1, 5)
                @test step_x(4, 2, mesh2d) == (2, 6)

                #Step y direction length more than 1
                @test step_y(x_size, 3, mesh2d) == (x_size - 3, 3)
                @test step_y(x_size - 3, 3, mesh2d) == (x_size - 6, x_size)
                @test step_y(4, 2, mesh2d) == (2, 6)
            end
        end
    end
end

@testset "3D AbstractMesh" begin
    @testset "Basic Functionality I" begin
        for x_size in 9:11
            for y_size in 9:11
                for z_size in 9:12
                    mesh3d = MetaMesh3D(x_size, y_size, z_size)

                    @test linear_size(mesh3d) == x_size * y_size * z_size

                    #Step x direction
                    @test step_x(1, mesh3d) == (y_size, 2)
                    @test step_x(y_size, mesh3d) == (y_size - 1, 1)
                    @test step_x(3, mesh3d) == (2, 4)

                    #Step y direction
                    @test step_y(1, mesh3d) == (x_size, 2)
                    @test step_y(x_size, mesh3d) == (x_size - 1, 1)
                    @test step_y(3, mesh3d) == (2, 4)

                    #Step z direction
                    @test step_z(1, mesh3d) == (z_size, 2)
                    @test step_z(z_size, mesh3d) == (z_size - 1, 1)
                    @test step_z(3, mesh3d) == (2, 4)

                    #Linear indexing back and forth
                    @test linear_indexing(3, 1, 1, mesh3d) == 3
                    @test linear_indexing(3, 2, 1, mesh3d) == x_size + 3

                    @test linear_indexing(x_size * y_size, mesh3d) == (x_size, y_size, 1)
                    @test linear_indexing(x_size, y_size, 1, mesh3d) == x_size * y_size

                    @test linear_indexing(x_size + 3, mesh3d) == (3, 2, 1)
                    @test linear_indexing(2 * x_size + 5, mesh3d) == (5, 3, 1)

                    @test linear_indexing(x_size * y_size * z_size, mesh3d) ==
                          (x_size, y_size, z_size)
                    @test linear_indexing(x_size * y_size * (z_size - 2), mesh3d) ==
                          (x_size, y_size, z_size - 2)

                    @test linear_indexing(x_size, y_size, 4, mesh3d) == x_size * y_size * 4
                    @test linear_indexing(x_size, 3, z_size, mesh3d) ==
                          x_size * y_size * (z_size - 1) + 3 * x_size

                    @test linear_indexing(x_size * y_size * 4, mesh3d) ==
                          (x_size, y_size, 4)
                    @test linear_indexing(x_size * y_size * (z_size - 1) + 3 * x_size,
                        mesh3d) == (x_size, 3, z_size)

                    #Step x direction length more than 1
                    @test step_x(1, 2, mesh3d) == (y_size - 1, 3)
                    @test step_x(3, 2, mesh3d) == (1, 5)
                    @test step_x(4, 2, mesh3d) == (2, 6)

                    #Step y direction length more than 1
                    @test step_y(x_size, 3, mesh3d) == (x_size - 3, 3)
                    @test step_y(x_size - 3, 3, mesh3d) == (x_size - 6, x_size)
                    @test step_y(4, 2, mesh3d) == (2, 6)

                    #Step z direction length more than 1
                    @test step_z(z_size, 3, mesh3d) == (z_size - 3, 3)
                    @test step_z(z_size - 3, 3, mesh3d) == (z_size - 6, z_size)
                    @test step_z(4, 2, mesh3d) == (2, 6)
                end
            end
        end
    end
end