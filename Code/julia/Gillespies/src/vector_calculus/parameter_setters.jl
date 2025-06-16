
"""
A function to take the Dict of parameter assignments, and place it into a
`@parameters` from ModelingToolkit

Using some of the metaprograming features found in Julia, we basically construct
the line of code that would say `@parameters x=1 y=2 ...` at parse time,
so it can be compiled into the line of code it would otherwise be.

Additionally we replace the existing dictionary that was fed into it.
"""
macro parameterification(expr)
    name = string(expr)
    out = quote dict_to_parameters($expr, $(Symbol(name))) end
    return esc(out)
end


"""
The function called by the macro,
I lost literal hours not noticing that I was missing a '...'
"""
function dict_to_parameters(D :: Dict{<: Any, <: Number}, name)
    syms  = [ :($k=$v) for (k, v) ∈ D ]
    expr₁ = quote @parameters $(syms...) end

    new_dict = [Meta.parse(replace("$val", "=" => "=>")) for val ∈ syms]
    expr₂ = quote new_dict = Dict(eval.($new_dict)) end

    expr₃ = quote push!($(name), new_dict...) end

    return quote $expr₁; $expr₂; $expr₃ end
end
