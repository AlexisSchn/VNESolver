# Automatically load the library globally when the module is used
function __init__()
    Libdl.dlopen("libkahip.so")
end



function kaffpa_core(
        n::Int32, 
        xadj::Vector{Int32}, 
        adjncy::Vector{Int32},
        nparts::Int32, 
        imbalance::Float64, 
        suppress_output::Bool, 
        seed::Int32, 
        mode::Int32
    )
    
    edgecut = Ref{Int32}(0)
    part = Vector{Int32}(undef, n)

    ccall(
        (:kaffpa, "libkahip.so"), 
        Cvoid,  
        ( 
            Ptr{Cint}, 
            Ptr{Cint}, 
            Ptr{Cint}, 
            Ptr{Cint}, 
            Ptr{Cint}, 
            Ptr{Cint}, 
            Ptr{Cdouble}, 
            Bool, 
            Cint, 
            Cint, 
            Ptr{Cint}, 
            Ptr{Cint}
        ),
        Ref(n), 
        C_NULL, 
        xadj, 
        C_NULL, 
        adjncy, 
        Ref(nparts), 
        Ref(imbalance), 
        suppress_output, 
        seed, 
        mode, 
        edgecut, 
        part
    )

    return Int(edgecut[]), part .+ 1 
end


function partition_kahip(g::Graph; nb_clusters=3, imbalance=0.1, seed=0, mode=0)
    println("Using KaHIP to partition the network. Warning: KaHIP in 3.25 necessary.")
    xadj, adjncy = graph_to_csr(g)
    n = Int32(nv(g))
    _, part = kaffpa_core(n, xadj, adjncy, Int32(nb_clusters), Float64(imbalance), true, Int32(seed), Int32(mode))
    return part
end