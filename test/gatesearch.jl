using GadgetSearcher, Test
using GadgetSearcher: match_gate

@testset "match gate" begin
    @logic NOR(x, y, z) := (z == ~(x | y))
    @logic XOR(x, y, z) := (z == x ⊻ y)
    graph_nor = simplegraph([(1, 2), (2, 3), (3, 4), (4, 5), (5, 1), (2, 4)])
    pins = [1, 3, 2]
    weights = ones(Int, 5)
    pw = PWGraph(graph_nor, pins, weights)
    @test GadgetSearcher.implements_gate(pw, NOR)
    @test !GadgetSearcher.implements_gate(pw, XOR)
end