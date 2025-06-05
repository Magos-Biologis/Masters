module Gillespies

# greet() = print("Hello World!")
# export StochasticSimuluationAlgorithm

using LinearAlgebra
using Random, Distributions

export ParamatersSSA,
       NovelSSA,
       ODESSA,
       IterationOfSSA,
       StepIterater


include("parameters.jl")
include("ssa.jl")

end # module Gillespies
