using SpecTools
using Test, Random

@testset "Base types" begin
  Random.seed!(2151)
  transitions = SpecTools.make_linear_levels(30) |> SpecTools.make_linear_transitions
  @test match_encoding(transitions, [0, 0, 1, 0]) == 1
  # test the matrix-ification of transitions
  tf = features(transitions)
  @test tf[:,1] == [200.0, 1.0, 0.0, 0.0, 1.0, 0.0] 
  # test grabbing the upper state energy levels
  @test upper_state_qnos(transitions[1]) == [1, 0]
end

@testset "Graph types" begin
  levels = SpecTools.make_linear_levels(20)
  transitions = SpecTools.make_linear_transitions(levels)
  # test constructing the graph objects
  proj_graph = ProjectedSG(levels, transitions)
  @test proj_graph isa SpecTools.ProjectedSG
  bi_graph = BipartiteSG(levels, transitions)
  @test bi_graph isa SpecTools.BipartiteSG
end

@testset "Lineshapes" begin
  @test Gaussian(1., 0., 0.5)(0.5) ≈ 0.48394144903
end

@testset "Partition structs" begin
  interp_obj = read_qpart("test.qpart")
  @test interp_obj isa SpecTools.InterpPartitionFunction
end

@testset "Partition function computation" begin
  @test partition_function(10_000, 1, 300f0) == 0.9984015312231446
end

