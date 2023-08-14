export gridsearch, gridsearch3, gridsearch4

function gridsearch4(target_matrix::AbstractArray{<:Tropical,4})
    locs0 = [(1.0, 3.0), (3.0, 1.0), (3.0, 5.0), (5.0, 3.0)]
    for d1=0:1, d2=0:1, d3=0:1, d4=0:1
        M1 = GenericTensorNetworks.content.(target_matrix)
        for i=0:(1<<9)-1
            locs = [locs0[1] .+ (0, d1), locs0[2] .+ (d2, 0), locs0[3] .+ (d3, 0), locs0[4] .+ (0, d4)]
            d1 == 1 && push!(locs, (0.0, 0.0))
            d2 == 1 && push!(locs, (5.0, 0.0))
            d3 == 1 && push!(locs, (5.0, 5.0))
            d4 == 1 && push!(locs, (0.0, 5.0))
            for k=0:8
                if (i>>k) & 1 == 1
                    x, y = k%3+2.0, k÷3+2.0
                    @assert x>=2 && x<=4
                    @assert y>=2 && y<=4
                    push!(locs, (x, y))
                end
            end
            g = unitdisk_graph(locs, 1.6)
            M2 = GenericTensorNetworks.content.(mis_compact_tropical_tensor(g, [1,2,3,4]))
            if is_diff_by_const(M1, M2)[1]
                @show locs
            end
        end
    end
    return nothing
end

function gridsearch3(target_matrix::AbstractArray{<:Tropical,3})
    M1 = GenericTensorNetworks.content.(target_matrix)
    locs0 = [(1.0, 3.0), (3.0, 1.0), (3.0, 5.0)]
    for i=0:(1<<9)-1
        locs = copy(locs0)
        for k=0:8
            if (i>>k) & 1 == 1
                x, y = k%3+2.0, k÷3+2.0
                @assert x>=2 && x<=4
                @assert y>=2 && y<=4
                push!(locs, (x, y))
            end
        end
        g = unitdisk_graph(locs, 1.6)
        M2 = GenericTensorNetworks.content.(mis_compact_tropical_tensor(g, [1,2,3]))
        if is_diff_by_const(M1, M2)[1]
            @show locs
        end
    end
    return nothing
end

function gridsearch(match, nrow, ncol, boundary_locs::T) where T<:Vector
    N = length(boundary_locs)
    res = T[]
    # check boundary locs
    for i=0:(1<<((nrow-2)*(ncol-2)))-1
        locs = copy(boundary_locs)
        for k=0:(nrow-2)*(ncol-2)-1
            if (i>>k) & 1 == 1
                y, x = k%(ncol-2)+2.0, k÷(ncol-2)+2.0
                @assert x>=2 && x<=nrow-1
                @assert y>=2 && y<=ncol-1
                push!(locs, (x, y))
            end
        end
        g = unitdisk_graph(locs, 1.6)
        M2 = GenericTensorNetworks.content.(mis_compact_tropical_tensor(g, collect(1:N)))
        if match(M2)
            @show locs
            push!(res, locs)
        end
    end
    return res
end

#        bi       si
#        ●        ●
#  ai x              ● ai
#
#  co ●              ● ci
#        x        ●
#        bi       so
function gate_multiplier()
    m = zeros(Bool, 2, 2, 2, 2, 2, 2)
    for ai=0:1, bi=0:1, ci=0:1, si=0:1
        lhs = (ai & bi) + si + ci
        so = lhs % 2
        co = lhs ÷ 2
        m[bi+1,si+1,ai+1,co+1,ci+1,so+1] = true
    end
    return m
end

export gridsearch_multiplier
function gridsearch_multiplier(nrow, ncol, x1, x2, y1, y2)
    mm = gate_multiplier()
    @show mm
    gridsearch(nrow, ncol, [(1, y1), (1, y2), (x1, ncol), (x2, 1), (x2, ncol), (nrow, y2)]) do m
        for i=1:length(m)
            !(mm[i]) == isinf(m[i]) || return false
        end
        return true
    end
end

export gridsearch_basicgate
function gate_basic(op)
    m = zeros(Bool, 2, 2, 2)
    for ai=0:1, bi=0:1
        ci = if op == :and
            ai & bi
        elseif op == :or
            ai | bi
        else
            @assert op == :xor
            ai ⊻ bi
        end
        m[ai+1, bi+1, ci+1] = true
    end
    return m
end

function gridsearch_basicgate(nrow, ncol, gatename)
    mm = gate_basic(gatename)
    @show gatename
    for locations in boundary_locations((nrow, ncol), 3)
        @show locations
        res = gridsearch(nrow, ncol, locations) do m
            for i=1:length(m)
                !(mm[i]) == isinf(m[i]) || return false
            end
            @show m
            return true
        end
        if !isempty(res)
            return res
        end
    end
end

export gridsearch_fulladder
function gate_fulladder()
    m = zeros(Bool, 2, 2, 2, 2, 2)
    for ai=0:1, bi=0:1, ci=0:1
        lhs = ai + bi + ci
        so = lhs % 2
        co = lhs ÷ 2
        m[ai+1, bi+1, ci+1, so+1, co+1] = true
    end
    return m
end

function gridsearch_fulladder(nrow, ncol)
    mm = gate_fulladder()
    #[(1, y1), (x1, ncol), (x2, ncol), (nrow, y1), (x2, 1)]
    for locations in boundary_locations((nrow, ncol), 5)
        @show locations
        res = gridsearch(nrow, ncol, locations) do m
            for i=1:length(m)
                !(mm[i]) == isinf(m[i]) || return false
            end
            return true
        end
        if !isempty(res)
            return res
        end
    end
end

function boundary_locations(sz, k)
    pool = Tuple{Int,Int}[]
    for i=1:sz[1]
        push!(pool, (i,1))
        push!(pool, (i,sz[2]))
    end
    for j=2:sz[2]-1
        push!(pool, (1,j))
        push!(pool, (sz[1],j))
    end
    iter_boundarys!(k, Tuple{Int,Int}[], pool, Vector{Tuple{Int,Int}}[])
end

function iter_boundarys!(n::Int, blist, pool, res)
    if n == 0
        push!(res, blist)
    end
    for loc in pool
        rem = filter(x->distance(loc, x) > 1.5, pool)
        iter_boundarys!(n-1, [blist..., loc], rem, res)
    end
    return res
end
distance(x, y) = sqrt(sum(abs2, x .- y))

export weighted_gridsearch
function wmissize(gp, weights)
    contractf(x->GenericTensorNetworks.TropicalF64(weights[x[1]]), gp)
end
function weighted_gridsearch()
    g0 = simplegraph([(1,5), (5,6), (6,7), (7,4), (2,8), (8,9), (9,10), (10, 3)])
    gp0 = Independence(g0, openvertices=[1,2,3,4])
    M1 = wmissize(gp0, vcat(ones(4), fill(2.0, nv(g0)-4)))
    locs0 = [(1.0, 3.0), (3.0, 1.0), (3.0, 5.0), (5.0, 3.0)]
    for d1=0:1, d2=0:1, d3=0:1, d4=0:1
        for i=0:(1<<9)-1
            locs = [locs0[1] .+ (0, d1), locs0[2] .+ (d2, 0), locs0[3] .+ (d3, 0), locs0[4] .+ (0, d4)]
            d1 == 1 && push!(locs, (0.0, 0.0))
            d2 == 1 && push!(locs, (5.0, 0.0))
            d3 == 1 && push!(locs, (5.0, 5.0))
            d4 == 1 && push!(locs, (0.0, 5.0))
            for k=0:8
                if (i>>k) & 1 == 1
                    x, y = k%3+2.0, k÷3+2.0
                    @assert x>=2 && x<=4
                    @assert y>=2 && y<=4
                    push!(locs, (x, y))
                end
            end
            g = unitdisk_graph(locs, 1.6)
            gp = Independence(g, openvertices=[1,2,3,4])
            M2 = wmissize(gp, vcat(ones(4), fill(2.0, length(locs)-4)))
            if all((maximum(M1) .== M1) .== (M2 .== maximum(M2)))
                @show locs
            end
        end
    end
    return nothing
end

function weighted_gridsearch3()
    g0 = simplegraph([(1,5), (5,6), (6,7), (7,4), (2,8), (8,9), (9,10), (10, 3)])
    gp0 = Independence(g0, openvertices=[1,2,3,4])
    M1 = wmissize(gp0, vcat(ones(4), fill(2.0, nv(g0)-4)))
    locs0 = [(1.0, 3.0), (3.0, 1.0), (3.0, 5.0), (5.0, 3.0)]
    for d1=0:1, d2=0:1, d3=0:1, d4=0:1
        for i=0:(1<<9)-1
            locs = [locs0[1] .+ (0, d1), locs0[2] .+ (d2, 0), locs0[3] .+ (d3, 0), locs0[4] .+ (0, d4)]
            d1 == 1 && push!(locs, (0.0, 0.0))
            d2 == 1 && push!(locs, (5.0, 0.0))
            d3 == 1 && push!(locs, (5.0, 5.0))
            d4 == 1 && push!(locs, (0.0, 5.0))
            for k=0:8
                if (i>>k) & 1 == 1
                    x, y = k%3+2.0, k÷3+2.0
                    @assert x>=2 && x<=4
                    @assert y>=2 && y<=4
                    push!(locs, (x, y))
                end
            end
            g = unitdisk_graph(locs, 1.6)
            gp = Independence(g, openvertices=[1,2,3,4])
            M2 = wmissize(gp, vcat(ones(4), fill(2.0, length(locs)-4)))
            if all((maximum(M1) .== M1) .== (M2 .== maximum(M2)))
                @show locs
            end
        end
    end
    return nothing
end

# the first N-1 dimensions are inputs
function get_truth_table(t::AbstractArray{T,N}) where {T,N}
    nin = 2^(N-1)
    m = reshape(t, nin, 2)
    res = zeros(Int, nin)
    MAX = maximum(m)
    for i=1:nin
        outval = findall(==(MAX), view(m, i, :))
        if length(outval) != 1
            return false, res
        end
        res[i] = outval[]-1
    end
    return true, res
end

function detect_logic_gate(g::Graph; inputs, output, weights=NoWeight())
    ixs = graph2ixs(g)
    size_dict = Dict([i=>2 for i=1:nv(g)])
    m = mis_compact_tropical_tensor(ixs, [inputs..., output], size_dict; compact=false, weights)
    #m = solve(IndependentSet(g; openvertices=[inputs..., output], weights), SizeMax())
    sig, table = get_truth_table(m)
    return sig, table
end

function int2graph(nv::Int, b::Int)
    g = SimpleGraph(nv)
    k = 0
    for i=1:nv
        for j=i+1:nv
            k += 1
            if !iszero(b>>(k-1) & 1)
                add_edge!(g, i, j)
            end
        end
    end
    return g
end

export search_weighted
function search_weighted(nv::Int; D::Int=1)
    innergraphfile = project_relative_path("data", "graph$nv.g6")
    graphs = load_g6(innergraphfile)
    for iw = 0:D^nv-1
        weights = [(iw÷(D^i) % D)+1 for i=0:nv-1]
        @show iw
        for g in graphs
            if is_connected(g)
                for inputs in Combinatorics.combinations(1:nv, 2)
                    for output in setdiff(1:nv, inputs)
                        sig, table = detect_logic_gate(g; inputs, output, weights)
                        if sig
                            #if table ∉ [[0, 0, 0, 0], [1,1,1,1], [1,1,0,0], [0,0,1,1],
                                    #[1, 0, 0, 0], [0,1,0,0], [0,0,1,0], [1,0,1,0], [0,1,0,1]]
                            if table == [0, 0, 0, 1]
                                nin = length(inputs)
                                ts = join([bitstring(i-1)[end-nin+1:end]*" => $(v)" for (i, v) in enumerate(table)], "\n")
                                @info "find logic gate, truth table is:\n$(ts)"
                                println([edges(g)...])
                                println(inputs, output)
                                println(weights)
                            end
                        end
                    end
                end
            end
        end
    end
end

# c = a ⊻ b