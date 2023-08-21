using Yao
using Yao.ConstGate: P1
using LinearAlgebra
using Graphs

function RydbergH(graph::SimpleGraph, U::Real, Δs::AbstractVector{<:Real}, Ωs::AbstractVector)
    Ωx = real.(Ωs)
    Ωy = imag.(Ωs)
    n = nv(graph)
    h = sum([U * kron(n, e.src=>P1, e.dst=>P1) for e in edges(graph)]) -
        sum([Δ * put(n, v=>P1) for (Δ, v) in zip(Δs, vertices(graph))]) +
        sum([Ω * put(n, v=>X) for (Ω, v) in zip(Ωx, vertices(graph))]) + 
        sum([Ω * put(n, v=>Y) for (Ω, v) in zip(Ωy, vertices(graph))])
    return h
end

function g5()
    g = SimpleGraph(5)
    for (i, j) in [(1, 2), (2, 3), (3, 4), (4, 5), (5, 1), (2, 5)]
        add_edge!(g, i, j)
    end
    return g
end

function effectiveH(z, H::AbstractMatrix, subspace::AbstractMatrix{T}) where T
    P = subspace * subspace'
    @assert P^2 ≈ P
    Q = I - P
    PHQ = P * H * Q
    return P * H * P + PHQ * pinv(z * Q - Q * H * Q) * PHQ'
end

@testset "effectiveH" begin
    n = 5
    h = EasyBuild.heisenberg(n; periodic=false)
    In = Matrix{ComplexF64}(I, 1<<n, 1<<n)
    E1 = eigvals(Matrix(h))
    PSI = In[:, 1:5]
    heff = effectiveH(E1[1], mat(h), PSI)
    E2 = eigvals(heff)
    heff3 = PSI' * heff * PSI
    E3 = eigvals(heff3)
    @test E1[1] ≈ E2[1]
    @test E1[1] ≈ E3[1]
end

# add up two graphs
function glue_graphs(g::SimpleGraph, wg::AbstractVector{T}, h::SimpleGraph, wh::AbstractVector{T}, gluepoints::Dict) where T
    # copy g
    graph = SimpleGraph(nv(g)+nv(h)-length(gluepoints))
    for e in edges(g)
        add_edge!(graph, e)
    end
    weights = zeros(T, nv(graph))
    weights[1:nv(g)] .= wg

    # add h
    vmap = Dict{Int,Int}()
    for (i, j) in gluepoints
        vmap[j] = i
        weights[i] += wh[j]
    end
    for (j, i) in zip(setdiff(vertices(h), values(gluepoints)), nv(g)+1:nv(graph))
        vmap[j] = i
        weights[i] = wh[j]
    end
    for e in edges(h)
        add_edge!(graph, vmap[e.src], vmap[e.dst])
    end

    return graph, weights
end

@testset "glue graphs" begin
    g = g5()
    h = g5()
    graph, weights = glue_graphs(g, ones(Int, 5), h, ones(Int, 5), Dict(3=>1))
    h = RydbergH(graph, 1000, weights, zeros(Int, nv(graph)))
    E = eigvals(Matrix(h))
    @test all(i->E[i] ≈ -4, 1:8)
    @test !(E[9] ≈ -4)
end

graph = g5()
Δs = 1e-2 * randn(5) .+ 1
Ωs = 1e-2 * randn(5)
h0 = RydbergH(graph, 1000.0, ones(Float64, 5), zeros(Float64, 5))
eig = eigen(Hermitian(Matrix(h0)))
# the lowest 4 energy levels
subspace = eig.vectors[:, 1:4]
h = RydbergH(graph, 1000.0, Δs, Ωs)
heff = subspace' * effectiveH(eig.values[1], mat(h), subspace) * subspace

Δs2 = 1e-2 * randn(5) .+ 1
Ωs2 = 1e-2 * randn(5)
# the lowest 4 energy levels
h2 = RydbergH(graph, 1000.0, Δs2, Ωs2)
heff2 = subspace' * effectiveH(eig.values[1], mat(h2), subspace) * subspace

# add up two effective Hamiltonian
effective123 = Matrix(put(3, (1, 2)=>matblock(heff)) + put(3, (2, 3)=>matblock(heff2)))

using Yao.YaoArrayRegister.StaticArrays
# creating mask
ggraph0, gweights0 = glue_graphs(graph, [SVector(a, b) for (a, b) in zip(zeros(Int, nv(graph)), ones(Int, nv(graph)))],
                                graph, [SVector(a, b) for (a, b) in zip(zeros(Int, nv(graph)), ones(Int, nv(graph)))],
                                Dict(3=>1))
gh0 = RydbergH(ggraph, 1000.0, getindex.(gweights0, 2), getindex.(gweights0, 1))
geig = eigen(Hermitian(Matrix(gh0)))
gsubspace = geig.vectors[:, 1:8]

ggraph, gweights = glue_graphs(graph, [SVector(a, b) for (a, b) in zip(Ωs, Δs)],
                                graph, [SVector(a, b) for (a, b) in zip(Ωs2, Δs2)],
                                Dict(3=>1))
h3 = RydbergH(ggraph, 1000.0, getindex.(gweights, 2), getindex.(gweights, 1))
effectiveabc = gsubspace' * effectiveH(geig.values[1], mat(h3), gsubspace) * gsubspace

effective123 .- effectiveabc
effective123