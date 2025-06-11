
# As this whole thing becomes a lot easier to manage with the structs concretely
# defined, I'm gonna go overboard
"""
This struct is for relevant counts gained from the `model_counts` function
"""
struct ParVarStruct{T<:Integer}
    par::T
    var::T
end

"""
This function just returns some integers from the name of the model,
for the sake of the parameter and variable setting.
"""
function ModelCounts(model::AbstractString) :: ParVarStruct #, is_ode::Bool)
    par_size = match(r"(\d+)Par", model).captures[1]
    var_size = match(r"(\d+)Var", model).captures[1]
    output = ParVarStruct( parse(Int, par_size), parse(Int, var_size) )
    return output
end


