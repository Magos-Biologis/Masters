module Gillespies

using LinearAlgebra

using Random
using Distributions



export ParametersSSA,
       NovelSSA,
       ODESSA,
       load_propensity_stuff,
       IterationOfSSA,
       StepIterater,
       vⱼ,
       propensities


# Figures out what the directory of this file is
const BASEDIR::String = @__DIR__
const PROPDIR::String = BASEDIR * "/propensities/"


include("parameters.jl")
include("loadprop.jl")

include("ssa.jl")





end # module Gillespies
