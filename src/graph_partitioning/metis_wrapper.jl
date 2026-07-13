

using METIS_jll

# Constants directly mapping to METIS documentation indices (0-indexed in C)
const METIS_NOPTIONS         = 40
const METIS_OPTION_CONTIG    = 11
const METIS_OPTION_SEED      = 12
const METIS_OPTION_UFACTOR   = 16

const METIS_OK               = Int32(1)



function metis_kway_core(n::Int32, xadj::Vector{Int32}, adjncy::Vector{Int32}, 
                         nparts::Int32, ufactor::Int32, seed::Int32, enforce_contiguous::Bool)
    
    ncon = Ref{Int32}(1)     # 1 balancing constraint
    objval = Ref{Int32}(0)   # Will store total edgecut value
    part = Vector{Int32}(undef, n)

    # Initialize METIS options vector
    options = Vector{Int32}(undef, METIS_NOPTIONS)
    ccall(
        (:METIS_SetDefaultOptions, METIS_jll.libmetis),
        Cvoid,
        (Ptr{Int32},),
        options
    )

    # Inject your requirements into the options vector (1-indexed offset for Julia)
    options[METIS_OPTION_SEED + 1] = seed
    options[METIS_OPTION_CONTIG + 1] = enforce_contiguous ? Int32(1) : Int32(0)
    
    # METIS expects ufactor as an integer representing max allowable imbalance percentage * 10
    # e.g., an imbalance tolerance of 3% is written as 30. Default is 30 (for kway).
    options[METIS_OPTION_UFACTOR + 1] = ufactor

    # Execute the partitioning graph execution
    status = ccall(
        (:METIS_PartGraphKway, METIS_jll.libmetis),
        Int32,
        (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
         Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
         Ptr{Float32}, Ptr{Float32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}),
        Ref(n), ncon, xadj, adjncy,
        C_NULL, C_NULL, C_NULL, # vwgt, vsize, adjwgt (NULL assumes uniform weight = 1)
        Ref(nparts), C_NULL, C_NULL, # tpwgts, ubvec (NULL assumes uniform balance targets)
        options, objval, part
    )

    if status != METIS_OK
        error("METIS partitioning failed with internal status code: $status")
    end

    return Int(objval[]), part .+ 1 # Vectorized adjustment back to 1-indexed Julia style
end


function partition_metis(g::Graph, nb_clusters::Int; imbalance_pct=3, seed=1, contiguous=true)
    if contiguous && !is_connected(g)
        error("METIS contiguous blocks feature requires the input graph to be fully connected!")
    end

    xadj, adjncy = graph_to_csr(g)
    n = Int32(nv(g))
    
    _, part = metis_kway_core(n, xadj, adjncy, Int32(nb_clusters), Int32(imbalance_pct * 10), Int32(seed), contiguous)
    return part
end