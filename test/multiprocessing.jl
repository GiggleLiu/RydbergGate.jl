using GadgetSearcher, Test
using GadgetSearcher: match_gate

@testset "multi processing" begin
    g5 = nonisomorphic_graphs(5)
    @logic NOR(x, y, z) := (z == ~(x | y))
    pws = Dict(["$name-pins431"=>PWGraph(g, [4, 3, 1], ones(Int, 5)) for (name, g) in g5])
    dirname = mktempdir()
    search_gate(pws, NOR; allow_degeneracy=false, multiprocess=true, output_folder=dirname)
    filename = joinpath(dirname, readdir(dirname)[1])
    pwgraph = load_pwgraph(filename)
    expected = PWGraph(simplegraph([(1, 3), (1, 4), (1, 5), (2, 4), (2, 5), (3, 5)]), [4, 3, 1], [1, 1, 1, 1, 1])
    @test pwgraph.graph == expected.graph
    @test pwgraph.weights == expected.weights
    @test pwgraph.pins == expected.pins
    @test pwgraph == expected
end