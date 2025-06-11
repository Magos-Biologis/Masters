"""
A Julia-fied version of my python implementation of the Gillespie Stochastic
Simulation Algorithm.
"""
function SSA(aⱼ:: Vector{T})::SSAOutputStruct where T<:Real
    a₀::Real = sum( aⱼ )
    output = SSAOutputStruct(0, 0.)
    if a₀ <= 0.
        println("womp womp")
        return output
    end

    r = rand(Uniform(eps(Float64),1), 2)
    rₐ₀::Float64 = r[2] * a₀
    output.τ = log(1 / r[1]) / a₀

    for sⱼ ∈  cumsum( aⱼ)
        output.j += 1
        sⱼ >= rₐ₀ && break ## The && is a 'short circuiting boolean AND'
    end

    return output
end



