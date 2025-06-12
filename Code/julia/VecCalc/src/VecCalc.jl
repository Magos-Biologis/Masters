module VecCalc

const BASEDIR::AbstractString = @__DIR__


using Symbolics #-> Just in case, though it may be reimported by ModelingToolkit
import Symbolics: jacobian

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D



include("common.jl")

include("odes/SimpleChemical1Par2Var.jl")

end
