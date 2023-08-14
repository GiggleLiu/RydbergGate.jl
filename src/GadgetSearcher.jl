module GadgetSearcher

using GenericTensorNetworks
using GenericTensorNetworks.Graphs, GenericTensorNetworks.OMEinsum
using DocStringExtensions
using Distributed
using MLStyle: @match
import GraphIO, JSON3

export nonisomorphic_graphs
export multiprocess_run
export search_gate
export @logic
export PWGraph, simplegraph


include("utils.jl")
include("gridgraph.jl")
include("gatesearch.jl")
include("multiprocessing.jl")

end
