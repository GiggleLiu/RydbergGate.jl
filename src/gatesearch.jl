"""
$TYPEDSIGNATURES

Generate the tensor for boolean logic.

### Example
```jldoctest setup=:(using GadgetSearcher)
julia> @logic T(x, y) := x < y
2×2 BitMatrix:
 0  1
 0  0

julia> @logic T(x, y) := x ⊻ y
2×2 BitMatrix:
 0  1
 1  0

julia> @logic T(x, y) := x & y
2×2 BitMatrix:
 0  0
 0  1
```

The returned value is directly assigned to `T`.
"""
macro logic(ex)
    @match ex begin
        :($name($(vars...)) := $expr) => begin
            nvar = length(vars)
            Expr(:block, [:($s = reshape([false, true], $(([k == i ? 2 : 1 for k in 1:nvar]...,)))) for (i, s) in enumerate(vars)]...,
                :($(esc(name)) = @. $expr)
            )
        end
    end
end

struct PWGraph{WT<:Union{NoWeight, Vector{Int}}}
    graph::SimpleGraph{Int}
    pins::Vector{Int}
    weights::WT
end

# function add_pins_and_weights(name::String, graph::SimpleGraph, weight_range)
#     d = Dict{String, PWGraph}()
#     for (k, weight_index) in enumerate(CartesianIndices(ntuple(i->length(weight_range), nv(graph))))
#         weights = [weight_range[k] for k in weight_index.I]
#         d[name * "-w$k"] = PWGraph(graph, pins, weights)
#     end
#     return d
# end

function search_gate(graphs::Dict{String, PWGraph{WT}}, target_gate; outputfolder, overwrite=false, multiprocess=false) where WT
    # check if the existing computing result exists
    if outputfolder !== nothing && isdir(outputfile) && !overwrite
        error("Non-empty folder without overwriting: $outputfolder")
    end
    if multiprocess
        # parallel run
        @info "The number of uncomputed inner graphs = $(length(uncomputed)), number of processes = $(nprocs())"
        multiprocess_run(graphs) do job
            name, graph = job
            _check_and_write(name, pw, target_gate)
        end
        return nothing
    else
        # serial run
        matched_graphs = PWGraph{WT}[]
        for (name, graph) in graphs
            @debug "Calculating graph: $name"
            _check_and_write(name, pw) && push!(matched_graphs, graph)
        end
        @info "# of matched graphs = $count"
        return matched_graphs
    end
end

# check if the `pw` implements the target gate, if so, write the result to the output folder.
function _check_and_write(name::String, pw::PWGraph, target_gate; outputfolder)
    res = implements_gate(pw, target_gate)
    if res && outputfolder !== nothing
        open(joinpath(outputfolder, "$name.dat"), "w") do f
            save_pwgraphs(f, pw)
        end
    end
    return res
end

# returns true if the `pw` graph implements the target gate.
function implements_gate(pw::PWGraph, target_gate::AbstractArray)
    n = nv(pw.graph)
    size_dict = Dict([i=>2 for i=1:n])
    ixs0 = graph2ixs(pw.graph)
    m = mis_compact_tropical_tensor(ixs0, pw.pins, size_dict; compact=false, weights=pw.weights)
    return match_gate(m, target_gate)
end

# test if the gate t matches the target gate tensor.
function match_gate(t::AbstractArray{<:Tropical, N}, gate_tensor::AbstractArray{Bool, N}) where N
    maxval = maximum(t)
    for i in eachindex(t)
        (t[i] == maxval) == gate_tensor[i] || return false
    end
    return true
end