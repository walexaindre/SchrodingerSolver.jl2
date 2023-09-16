export plot_inner
function plot_inner(value, mesh::SpaceTimeGrid2D,id::Int = 0)
    indexes = get_indexed_elements(mesh)
    indexes = [x[k] for x in indexes, k in 1:2]

    #absv = sum(abs2.(value), dims=2) |> Array
    absv = abs2.(value[:,1])

    fig = Figure()
    axe = Axis3(fig[1, 1])
    plot!(axe, indexes[:, 1], indexes[:, 2], absv[:, 1])
    save("fig$id.png",fig)
    display(fig)
end

function plot_inner(value, mesh::SpaceTimeGrid1D,id::Int = 0)
    indexes = get_indexed_elements(mesh)
    #indexes = [x[k] for x in indexes, k in 1:2]

    #absv = sum(abs2.(value), dims=2) |> Array
    absv = abs2.(value[:,1])

    fig = Figure()
    axe = Axis(fig[1, 1])

    @show typeof(absv) typeof(indexes)

    plot!(axe, indexes, absv)
    #save("fig$id.png",fig)
    display(fig)
end