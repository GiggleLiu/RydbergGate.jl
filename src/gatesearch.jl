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

struct PWGraph
    graph::SimpleGraph{Int}
    pins::Vector{Int}
    weights::Vector{Int}
end
Base.:(==)(a::PWGraph, b::PWGraph) = a.graph == b.graph && a.pins == b.pins && a.weights == b.weights

"""
$TYPEDSIGNATURES

Save a [`PWGraph`](@ref) object to a file.
"""
function save_pwgraph(filename::String, pw::PWGraph)
    open(filename, "w") do f
        JSON3.write(f, pw)
    end
end

"""
$TYPEDSIGNATURES

Load a [`PWGraph`](@ref) object from a file.
"""
function load_pwgraph(filename::String)
    return open(filename, "r") do f
        JSON3.read(f, PWGraph)
    end
end

# function add_pins_and_weights(name::String, graph::SimpleGraph, weight_range)
#     d = Dict{String, PWGraph}()
#     for (k, weight_index) in enumerate(CartesianIndices(ntuple(i->length(weight_range), nv(graph))))
#         weights = [weight_range[k] for k in weight_index.I]
#         d[name * "-w$k"] = PWGraph(graph, pins, weights)
#     end
#     return d
# end

function search_gate(graphs::Dict{String, PWGraph}, target_gate; output_folder=nothing, overwrite=false, multiprocess=false, allow_degeneracy=true)
    # check if the existing computing result exists
    if output_folder !== nothing && isdir(output_folder) && !isempty(readdir(output_folder)) && !overwrite
        error("Non-empty folder without overwriting: $output_folder")
    end
    if output_folder === nothing && multiprocess
        error("You should specify an output folder for a multi-processing run!")
    end
    if multiprocess
        # parallel run
        @info "The number of graphs = $(length(graphs)), number of processes = $(nprocs())"
        multiprocess_run(collect(graphs)) do job
            name, pw = job
            _check_and_write(name, pw, target_gate; allow_degeneracy, output_folder)
        end
        return nothing
    else
        # serial run
        matched_graphs = PWGraph[]
        for (name, pw) in graphs
            @debug "Calculating graph: $name"
            _check_and_write(name, pw, target_gate; allow_degeneracy, output_folder) && push!(matched_graphs, pw)
        end
        @info "# of matched graphs = $(length(matched_graphs))"
        return matched_graphs
    end
end

# check if the `pw` implements the target gate, if so, write the result to the output folder.
function _check_and_write(name::String, pw::PWGraph, target_gate; output_folder, allow_degeneracy)
    res = implements_gate(pw, target_gate; allow_degeneracy)
    if res && output_folder !== nothing
        save_pwgraph(joinpath(output_folder, "$name.dat"), pw)
    end
    return res
end

# returns true if the `pw` graph implements the target gate.
function implements_gate(pw::PWGraph, target_gate::AbstractArray; allow_degeneracy=true)
    n = nv(pw.graph)
    size_dict = Dict([i=>2 for i=1:n])
    ixs0 = graph2ixs(pw.graph)
    m = mis_tropical_tensor(ixs0, pw.pins, size_dict; compact=false, weights=pw.weights)
    res = match_gate(m, target_gate)
    if res && !allow_degeneracy
        # check degeneracy
        countings = solve(IndependentSet(pw.graph, optimizer=GreedyMethod(; nrepeat=1), openvertices=pw.pins; weights=pw.weights), CountingMax())
        maxval = maximum(x->x.n, countings)
        return all(x->x.n != maxval || x.c == 1, countings)
    else
        return res
    end
end

# test if the gate t matches the target gate tensor.
function match_gate(t::AbstractArray{<:Tropical, N}, gate_tensor::AbstractArray{Bool, N}) where N
    maxval = maximum(t)
    for i in eachindex(t)
        (t[i] == maxval) == gate_tensor[i] || return false
    end
    return true
end