using Distributed, Pkg
@everywhere using Distributed, Pkg
if nprocs() == 1
    addprocs(2)
end

using GadgetSearcher

@everywhere using GadgetSearcher
folder = GadgetSearcher.project_relative_path("data", ARGS[3])
if !isdir(folder)
    mkdir(folder)
end

counts = search_gate(parse(Int,ARGS[1]); cross_type=ARGS[3], outputfolder=folder, overwrite=parse(Bool,ARGS[2]), onlycount=false)
println("total count = $(sum(counts))")
