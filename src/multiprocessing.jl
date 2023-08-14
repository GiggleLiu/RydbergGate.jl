export multiprocess_run
export search_gadget_all_multiprocess

function do_work(f, jobs, results) # define work function everywhere
    while true
        job = take!(jobs)
        @info "running $job"
        res = f(job)
        put!(results, res)
    end
end

function multiprocess_run(func, inputs::AbstractVector{T}) where T
    n = length(inputs)
    jobs = RemoteChannel(()->Channel{T}(n));
    results = RemoteChannel(()->Channel{Any}(n));
    for i in 1:n
        put!(jobs, inputs[i])
    end
    for p in workers() # start tasks on the workers to process requests in parallel
        remote_do(do_work, p, func, jobs, results)
    end
    return Any[take!(results) for i=1:n]
end

function search_gadget_all_multiprocess(n; cross_type="simple", overwrite=false, outputfolder, onlycount=false, opt_connected=true, opt_inner_isometry=true, opt_can_cross=true)
    ngraph = nonisomorphicgraph_count[n]
    if !overwrite
        uncomputed = Int[]
        for i=1:ngraph
            f = joinpath(outputfolder, "results-$(n)$(i)-$(i).adj.dat")
            !isfile(f) && push!(uncomputed, i)
        end
    else
        uncomputed = collect(1:ngraph)
    end
    @info "The number of uncomputed inner graphs = $(length(uncomputed)), number of processes = $(nprocs())"
    #return
    innergraphfile = project_relative_path("data", "graph$n.g6")
    multiprocess_run(collect(zip(uncomputed, load_g6(innergraphfile)[uncomputed]))) do job
        i, graph = job
        res, count = search_gadget_and_write(graph, i; cross_type=cross_type, overwrite=overwrite, outputfolder=outputfolder,
        onlycount=onlycount, opt_connected=opt_connected, opt_inner_isometry=opt_inner_isometry, opt_can_cross=opt_can_cross)
        return count
    end
end