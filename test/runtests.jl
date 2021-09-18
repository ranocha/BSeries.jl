using Test
using BSeries

using Symbolics: Symbolics, @variables, Num


@testset "subsitution" begin
  @variables a1 a2 a31 a32
  a = OrderedDict{RootedTrees.RootedTree{Int, Vector{Int}}, Num}(
        rootedtree(Int[])     => 1,
        rootedtree([1])       => a1,
        rootedtree([1, 2])    => a2,
        rootedtree([1, 2, 2]) => a31,
        rootedtree([1, 2, 3]) => a32)

  @variables b1 b2 b31 b32
  b = OrderedDict{RootedTrees.RootedTree{Int, Vector{Int}}, Num}(
        rootedtree(Int[])     => 0,
        rootedtree([1])       => b1,
        rootedtree([1, 2])    => b2,
        rootedtree([1, 2, 2]) => b31,
        rootedtree([1, 2, 3]) => b32)

  t = rootedtree([1])
  @inferred substitute(b, a, t)
  @test isequal(substitute(b, a, t), a1 * b1)

  t = rootedtree([1, 2, 3])
  @inferred substitute(b, a, t)
  @test isequal(substitute(b, a, t), a1 * b32 + a32 * b1^3 + 2 * a2 * b1 * b2)
end
