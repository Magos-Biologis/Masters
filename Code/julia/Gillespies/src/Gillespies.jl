module Gillespies

using LinearAlgebra

using Random
using Distributions



export ParameterStructSSA,
       NovelSSA,
       ODESSA,
       load_propensity_stuff,
       IterationOfSSA,
       StepIterater,
       vⱼ,
       propensities


# Figures out what the directory of this file is
const BASEDIR::String = @__DIR__
const PROPDIR::String = joinpath(BASEDIR, "propensities")


include("common.jl")
include("loadprop.jl")

# include("engine.jl")
include("ssa.jl")



end # module Gillespies
