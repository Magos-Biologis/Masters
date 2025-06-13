"""
A Julia-fied version of my python implementation of the Gillespie Stochastic
Simulation Algorithm.

Made to function 'in-place' by just modifying the same preallocated array.
"""
function SSA!(output :: SSAOutput, aⱼ:: Vector{T}) :: SSAOutput where T<:Real
    a₀::Real = sum( aⱼ )
    output.j = 0;

    a₀ <= 0. ? ( println("womp womp");  return output) : nothing

    r = rand(Uniform(eps(Float64),1), 2)
    rₐ₀::Float64 = r[2] * a₀
    output.τ = log(1 / r[1]) / a₀

    for (i, sⱼ) ∈  enumerate(cumsum( aⱼ))
        ## The && is a 'short circuiting boolean AND'
        sⱼ >= rₐ₀ && (output.j = i; break)
        # sⱼ >= rₐ₀ ? output.j = i : nothing
    end

    return output
end



