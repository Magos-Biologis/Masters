import ModelingToolkit.Symbolics.variable as variable
import ModelingToolkit.Symbolics.variables as variables

"""
A function to take the Dict of parameter assignments, and place it into a
`@parameters` from ModelingToolkit

Using some of the metaprograming features found in Julia, we basically construct
the line of code that would say `@parameters x=1 y=2 ...` at parse time,
so it can be compiled into the line of code it would otherwise be.
"""
macro parameterification(dictionary)
    local dict = esc(dictionary)

    local list₁ = :( [ :($k=$v) for (k,v) ∈  $dictionary ] )
    local list₂ = :( variable.( keys($dict) ) .=> values($dict) )

    local expr₁ = :( $parameter_expansion($list₁) )
    local expr₂ = :( push!($dict, $list₂...) )

    return quote $(esc(:(eval($expr₁)))); $expr₂ end

end

function parameter_expansion(exprs :: V) where V <: Vector{Expr}
    expr = :( @parameters $(exprs...) )
    return quote $expr end
end



# if !(Meta.isexpr(dictionary, :dict) | typeof(dictionary.args[1]) <: Dict)
#     error("@parameterification must receive a literal Dict(...)")
# end

