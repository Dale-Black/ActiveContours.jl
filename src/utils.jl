function central_diff(u)
    du1 = diff(u; dims=1) / 2
    a = [du1[[1], :]; du1]
    a .+= [du1; du1[[end], :]]

    du2 = diff(u; dims=2) / 2
    b = [du2[[1], :]; du2]
    b .+= [du2; du2[[end], :]]

    return a, b
end
