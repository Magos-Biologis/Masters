module Gillespies

using LinearAlgebra

using Random
using Distributions


export ParameterTypeSSA,
       NovelStructSSA,
       DifferentialStructSSA,
       load_propensity_stuff,
       SSA,
       StepIterater
       # vⱼ,
       # propensities


# Figures out what the directory of this file is
const BASEDIR::AbstractString = @__DIR__
const PROPDIR::AbstractString = joinpath(BASEDIR, "propensities")


include("common.jl")
include("loadprop.jl")
include("ssa.jl")


# include("propensities/SimpleTwoChemicalSpecies.jl")


end # module Gillespies
