using Combinatorics

export search_gadget_gate, nand_mat, nor_mat, search_gadget_gate_all
export NOR, NAND, AbstractGate, XOR, AnyGate, AND

abstract type AbstractGate end
struct NOR <: AbstractGate end
struct NAND <: AbstractGate end
struct XOR <: AbstractGate end
struct AND <: AbstractGate end
struct AnyGate <: AbstractGate end

function search_gadget_gate(parentgraph::SimpleGraph, target_gate::AbstractGate; onlycount::Bool=false, opt_connected=true)
    res = SimpleGraph{Int}[]
    count = 0
    n = nv(parentgraph)
    size_dict = Dict([i=>2 for i=1:n])
    #isosets = isometric_sets(parentgraph)
    ixs0 = graph2ixs(parentgraph)
    symmetric12, symmetric23, symmetric13 = detect_symmetry(gatetensor(target_gate))
    @time for c=combinations(1:n, 3)
        (!opt_connected || is_connected(g)) || continue
        count += 1
        onlycount && continue
        m = mis_compact_tropical_tensor(ixs0, c, size_dict, false)
        M1 = GenericTensorNetworks.content.(m)
        if is_match_gate(M1, target_gate)[1] ||
            (!symmetric12 && is_match_gate(permutedims(M1,(2,1,3)), target_gate)) ||
            (!symmetric23 && is_match_gate(permutedims(M1,(1,3,2)), target_gate)) ||
            (!symmetric13 && is_match_gate(permutedims(M1,(3,2,1)), target_gate))
            @info "$c is good!"
            push!(res, g)
        end
    end
    return res, count
end

function is_match_gate(t, g::XOR)
    gmat = gatetensor(g)#content.(mis_compactify!(Tropical{Float64}.(gatetensor(g))))
    @show t
    for i=1:length(t)
        @show iszero(gmat[i]), isinf(t[i]), t[i]
        iszero(gmat[i]) == isinf(t[i]) || return false
    end
    return true
end

function is_match_gate(t, g::AND)
    gmat = gatetensor(g)#content.(mis_compactify!(Tropical{Float64}.(gatetensor(g))))
    for i=1:length(t)
        @show iszero(gmat[i]), isinf(t[i]), t[i]
        iszero(gmat[i]) == isinf(t[i]) || return false
    end
    return true
end

function is_match_gate(t, g)
    maxval = maximum(t)
    gmat = gatetensor(g)
    for i=1:length(t)
        (t[i] == maxval) == gmat[i] || return false
    end
    return true
end

function is_match_gate(t, g::AnyGate)
    maxval = maximum(t)
    gmat = gatetensor(g)
    for i=1:length(t)
        (t[i] == maxval) == gmat[i] || return false
    end
    return true
end

function gatetensor(::NOR)
    m = zeros(Bool, 2, 2, 2)
    for i=(false,true)
        for j=(false,true)
            m[i+1,j+1,!(i||j)+1] = 1
        end
    end
    return m
end

function gatetensor(::NAND)
    m = zeros(Bool, 2, 2, 2)
    for i=(false,true)
        for j=(false,true)
            m[i+1,j+1,!(i&&j)+1] = 1
        end
    end
    return m
end

function gatetensor(::XOR)
    m = zeros(Bool, 2, 2, 2)
    for i=(false,true)
        for j=(false,true)
            m[i+1,j+1,(i⊻j)+1] = 1
        end
    end
    return m
end

function gatetensor(::AND)
    m = zeros(Bool, 2, 2, 2)
    for i=(false,true)
        for j=(false,true)
            m[i+1,j+1,(i&&j)+1] = 1
        end
    end
    return m
end

function detect_symmetry(cmat::AbstractArray{T,3}) where T
    symmetric12 = permutedims(cmat, (2,1,3)) == cmat
    symmetric23 = permutedims(cmat, (1,3,2)) == cmat
    symmetric13 = permutedims(cmat, (3,2,1)) == cmat
    return symmetric12, symmetric23, symmetric13
end

function search_gadget_gate_all(n, target_gate; outputfolder, onlycount=false, overwrite=false, kwargs...)
    gatename = typeof(target_gate)
    if outputfolder!==nothing
        if isfile(joinpath(outputfolder, "results-$gatename-$(n).adj.dat")) && !overwrite
            println("Find existing computing results...")
            return
        end
    end
    if n==1
        res, count = search_gadget_gate(SimpleGraph(1), target_gate; onlycount=onlycount, kwargs...)
        return count
    end
    count = 0
    res = SimpleGraph{Int}[]
    iterateg6(project_relative_path("data", "graph$(n).g6")) do k, graph
        @info "Calculating n = $n, k=$k"
        resi, ci = search_gadget_gate(graph, target_gate; kwargs...)
        count += ci
        append!(res, resi)
    end
    if outputfolder !== nothing && !onlycount
        filename = joinpath(outputfolder, "results-$gatename-$n.adj.dat")
        _save_results(n, res, filename)
    end
    @info "# of configurations = $count"
    return res, count
end