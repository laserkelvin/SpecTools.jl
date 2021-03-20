using SpecTools
using Test

@testset "Energy Levels" begin
  # test the creation of the energy levels
  @test LinearLevel(109.2, 1, 1, 0) isa LinearLevel
  @test AsymTopLevel(5190256., 50, 20, 5, 3, 0) isa AsymTopLevel
end

@testset "Lineshapes" begin
  @test Gaussian(0., 1., 1.)(0.5) ≈ 0.14045374430
end
