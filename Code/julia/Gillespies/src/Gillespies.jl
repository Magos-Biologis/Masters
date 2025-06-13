__precompile__()
module Gillespies

using LinearAlgebra

using Random
import Distributions: Uniform

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

using StochasticDiffEq

using JSON
using NPZ


# The initial stuff is the structs of importance
export ParameterType,
       NovelStruct,
       DifferentialStruct,

       PropensityType,
       StepperStruct,

       # load_propensity_stuff,
       SSA!,
       StepIterator,
       SSAOutput,

       ### Vector Calculus stuff
       # SimpleChemical1Par2Var,


       LangevinEquation


       # SaveToNPZ


# Figures out what the directory of this file is
const BASEDIR::AbstractString = @__DIR__
const MODEDIR::AbstractString = joinpath(BASEDIR, "models")


include("common.jl")

include("stochastics/loadprop.jl")
include("stochastics/ssa.jl")
include("stochastics/stepper.jl")

include("vector_calculus/langevin.jl")

# The stochastic sim modules
# including it programatically cause why not
models = readdir(MODEDIR)
for model in models
    model_dir = joinpath(MODEDIR, model)
    include(joinpath(model_dir, "propensities.jl"))
    include(joinpath(model_dir, "ode.jl"))
end


end # module Gillespies
