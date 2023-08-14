project_relative_path(xs...) = normpath(joinpath(dirname(dirname(pathof(@__MODULE__))), xs...))

function ixs2graph(n::Int, ixs)
    g = SimpleGraph(n)
    for (i, j) in ixs
        add_edge!(g, i, j)
    end
    return g
end

function graph2ixs(g::SimpleGraph)
    return [[e.src, e.dst] for e in edges(g)]
end

# create a simple graph from edges
function mis_tropical_tensor(ixs::AbstractVector, iy::AbstractVector, size_dict; compact::Bool=true, weights=NoWeight())
    ixs = vcat([[i] for i in 1:length(size_dict)], ixs) # labels for vertex tensors
    optcode = optimize_fastgreedy(ixs, iy)  # specialized, faster approach
    xs = [length(ix)==1 ? GenericTensorNetworks.misv([Tropical(0.0), Tropical(weights isa NoWeight ? 0.0 : weights[ix[1]])]) : GenericTensorNetworks.misb(Tropical{Float64}, length(ix)) for ix in ixs]
    m = einsum(optcode, (xs...,), size_dict)
    compact && mis_compactify!(m)
    return m
end

# a faster greedy optimizer
function optimize_fastgreedy(ixs, iy)
    mask = trues(length(ixs))
    counts = Dict{Int, Int}()
    for ix in ixs
        for l in ix
            counts[l] = get(counts, l, 0) + 1
        end
    end
    for l in iy
        counts[l] = get(counts, l, 0) + 1
    end
    res = DynamicNestedEinsum{Int}(1)
    mask[1] = false
    for l in ixs[1]
        counts[l] -= 1
    end
    iout = ixs[1]
    for i=1:length(ixs)-1
        j = findfirst(1:length(ixs)) do j
            mask[j] && any(l->l ∈ iout, ixs[j])
        end
        if j === nothing
            j = findfirst(mask)
        end
        for l in ixs[j]
            counts[l] -= 1
        end
        iout_new = filter(l->counts[l]>0, iout ∪ ixs[j])
        res = DynamicNestedEinsum((res, DynamicNestedEinsum{Int}(j)), DynamicEinCode([iout, ixs[j]],i==length(ixs)-1 ? iy : iout_new))
        iout = iout_new
        mask[j] = false
    end
    return res
end


function mis_tropical_tensor(g::SimpleGraph, openvertices::AbstractVector, size_dict=Dict([i=>2 for i=1:nv(g)]); compact::Bool=true)
    ixs = graph2ixs(g)
    return mis_tropical_tensor(ixs, openvertices, size_dict; compact=compact)
end

# create a simple graph from edges
function simplegraph(edges)
    n = maximum(x->max(x...), edges)
    g = SimpleGraph(n)
    for (i,j) in edges
        add_edge!(g, i,j)
    end
    return g
end

# two tensors are different by a constant
function is_diff_by_const(t1::AbstractArray{T}, t2::AbstractArray{T}) where T
	x = NaN
	for (a, b) in zip(t1, t2)
		if isinf(a) && isinf(b)
			continue
		end
		if isinf(a) || isinf(b)
			return false, 0
		end
		if isnan(x)
			x = (a - b)
		elseif x != a - b
			return false, 0
		end
	end
	return true, x
end

"""
$TYPEDSIGNATURES

Enumerate all small size non-isomorphic graphs.
"""
function nonisomorphic_graphs(n::Int)
    if n==1
        return Dict("graph1"=>SimpleGraph(1))
    elseif n==0
        return Dict("graph1"=>SimpleGraph(0))
    elseif n > 8
        # please check: https://www.juliabloggers.com/the-return-of-the-graphs-and-an-interesting-puzzle/
        error("Data for non-isometric graphs with size n > 8 are not shiped in this package, please specify the `g6file` keyword argument.")
    end
    graphfile = project_relative_path("data", "graph$(n).g6")
    return open(graphfile) do f
        GraphIO.Graph6.loadgraph6_mult(f)
    end
end