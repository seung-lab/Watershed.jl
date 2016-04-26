VERSION >=v"0.4.0-dev+6521" && __precompile__()

module Watershed

# load dependencies
#include("steepestascent.jl")
# include("divideplateaus.jl")
# include("findbasins.jl")
# include("regiongraph.jl")
# include("mergeregions.jl")
# include("mst.jl")

# export steepestascent, divideplateaus!, findbasins, regiongraph, mergeregions, mst

include("watershed.jl")

end
