include("walls-and-chambers.jl")

function main()

    # Q, d = Quiver("1-2,1-3,1-4,2-3,3-4,"), [1, 1, 1, 1]
    Q, d = subspace_quiver(6), [1, 1, 1, 1, 1, 1, 2]


    @info "Quiver: $(Q)"
    @info "Dimension vector: $(d)"

    @info ""
    @info "Computing the VGIT fan. This might take a while..."
    fan = vgit_fan(Q, d; verbose=true)
    @info ""
    @info "Done. The f-vector of the fan is $(Oscar.f_vector(fan))."

    serialize("$(Q.adjacency),$(d)_fan.jls", fan)
end

main()
