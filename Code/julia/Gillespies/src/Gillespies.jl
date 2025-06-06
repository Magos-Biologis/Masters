module Gillespies

using LinearAlgebra

using Random
using Distributions

using JSON
using NPZ


export ParameterTypeSSA,
       NovelStructSSA,
       DifferentialStructSSA,

       # load_propensity_stuff,
       SSA,
       StepIterator,

       SaveToNPZ


# Figures out what the directory of this file is
const BASEDIR::AbstractString = @__DIR__
const PROPDIR::AbstractString = joinpath(BASEDIR, "propensities")


include("common.jl")

include("loadprop.jl")
include("ssa.jl")
include("stepper.jl")

include("saving.jl")


# The stochastic sim modules
# including it programatically cause why not
files = readdir(PROPDIR)
for file in files
    include(joinpath(PROPDIR, file))
end
# This seems to work a lot better than loading modules mid function,
# as loading the modules at runtime overwrites the module and/or
# prevents it from working
# Could be a Julia REPL thing, but if it works in there, it'll work
# elsewhere, so it's probably safer this way



end # module Gillespies
