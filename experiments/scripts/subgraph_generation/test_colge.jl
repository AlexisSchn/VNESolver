# to test the flow formulation

using CSV, DataFrames
using VNESolver

size = ARGS[1]
path_instances = "instances/$size/"
path_results = ARGS[2]

println("Path instances: $path_instances")
println("Path results: $path_results")

vns = read_virtuals_folder(path_instances*"vns/")
sns = read_substrates_folder(path_instances*"sns/")

# Warm up
dummy_vn = vns[1]
dummy_sn = sns[1]
dummy_instance = Instance(dummy_vn, dummy_sn)
_ = solve_subgraph_decomposition_better(dummy_instance)

results = AbstractSolverResult[]
 
# Loop
for vn in vns
    println("Doing vn $(vn.name)")
    for sn in sns
        println("   for sn $(sn.name)")

        instance = Instance(vn, sn)
        result = solve_subgraph_decomposition_better(instance)

        push!(results, result)        
        df_results = DataFrame(results)
        CSV.write("$path_results/results_$size.csv", df_results)
    end
end

println("Process complete. Saved results to $path_results/results_$size.csv")


