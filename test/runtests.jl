

using Test
using VNESolver

@testset "Flow formulation" begin
    instance = get_instance_from_folder("instances/1/")
    result = solve_flow_formulation(instance)

    @test result[:mapping_cost] == 108.
end