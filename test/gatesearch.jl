using GadgetSearcher, Test
using GadgetSearcher: match_gate

@testset "match gate" begin
    @logic NOR(x, y, z) := (z == ~(x | y))
    @logic XOR(x, y, z) := (z == x ⊻ y)
    graph_nor = simplegraph([(1, 2), (2, 3), (3, 4), (4, 5), (5, 1), (2, 4)])
    # this configuration contains degenerate space
    graph_nor2 = simplegraph([(1, 2), (2, 3), (3, 4), (4, 5), (5, 1)])
    pins = [1, 3, 2]
    weights = ones(Int, 5)
    pw = PWGraph(graph_nor, pins, weights)
    pw2 = PWGraph(graph_nor2, pins, weights)
    @test GadgetSearcher.implements_gate(pw, NOR; allow_degeneracy=false)
    @test !GadgetSearcher.implements_gate(pw2, NOR; allow_degeneracy=false)
    @test GadgetSearcher.implements_gate(pw2, NOR; allow_degeneracy=true)
    @test !GadgetSearcher.implements_gate(pw, XOR; allow_degeneracy=true)

    # check and write
    @test GadgetSearcher._check_and_write("", PWGraph(graph_nor, pins, ones(Int, 5)), NOR; output_folder=nothing, allow_degeneracy=false)

    # weighted graphs
    @logic OR(x, y, z) := (z == x | y)
    graph_or = simplegraph([(1, 2), (2, 3), (3, 4), (4, 5), (5, 1), (2, 4), (2, 6)])
    pins = [1, 3, 6]
    weights = ones(Int, 6)
    weights[2] = 2
    pw = PWGraph(graph_or, pins, weights)
    @test GadgetSearcher.implements_gate(pw, OR)
end

@testset "search gate" begin
    # loading all 5 vertex non-isomorphic graphs
    g5 = nonisomorphic_graphs(5)
    @logic NOR(x, y, z) := (z == ~(x | y))
    pws = Dict(["$name-pins123"=>PWGraph(g, [1, 2, 3], ones(Int, 5)) for (name, g) in g5])
    results = search_gate(pws, NOR; allow_degeneracy=false, multiprocess=false)
    @test length(results) == 0
    pws = Dict(["$name-pins431"=>PWGraph(g, [4, 3, 1], ones(Int, 5)) for (name, g) in g5])
    results = search_gate(pws, NOR; allow_degeneracy=false, multiprocess=false)
    @test length(results) == 1
    # saving and loading files
    pws = Dict(["$name-pins431"=>PWGraph(g, [4, 3, 1], ones(Int, 5)) for (name, g) in g5])
    dirname = mktempdir()
    results = search_gate(pws, NOR; allow_degeneracy=false, multiprocess=false, output_folder=dirname)
    @test length(results) == 1
    filename = joinpath(dirname, readdir(dirname)[1])
    pwgraph = load_pwgraph(filename)
    expected = PWGraph(simplegraph([(1, 3), (1, 4), (1, 5), (2, 4), (2, 5), (3, 5)]), [4, 3, 1], [1, 1, 1, 1, 1])
    @test pwgraph.graph == expected.graph
    @test pwgraph.weights == expected.weights
    @test pwgraph.pins == expected.pins
    @test pwgraph == expected
end