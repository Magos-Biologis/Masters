module Gillespies

# import Pkg

using LinearAlgebra

using Random
import Distributions: Uniform

using JSON
using NPZ


# The initial stuff is the structs of importance
export ParameterTypeSSA,
       NovelStructSSA,
       DifferentialStructSSA,

       PropensityType,
       StepperStruct,

       # load_propensity_stuff,
       SSA,
       StepIterator

       # SaveToNPZ


# Figures out what the directory of this file is
const BASEDIR::AbstractString = @__DIR__
const PROPDIR::AbstractString = joinpath(BASEDIR, "propensities")


include("common.jl")

include("loadprop.jl")
include("ssa.jl")
include("stepper.jl")

# The stochastic sim modules
# including it programatically cause why not
files = readdir(PROPDIR)
for file in files
    include(joinpath(PROPDIR, file))
    # Pkg.import(joinpath(PROPDIR, file))
end

# This seems to work a lot better than loading modules mid function,
# as loading the modules at runtime overwrites the module and/or
# prevents it from working
# Could be a Julia REPL thing, but if it works in there, it'll work
# elsewhere, so it's probably safer this way


# # We will just manually load them all actually,
# include("propensities/SimpleTwoChemicalSystem.jl")
#
# include("propensities/ODE3Par2Var.jl")
# include("propensities/ODE3Par3Var.jl")
#
# include("propensities/NovelFiveTwo.jl")
# include("propensities/NovelFiveThree.jl")


end # module Gillespies
