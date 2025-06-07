"""
A Julia-fied version of my python implementation of the Gillespie Stochastic
Simulation Algorithm.
"""
function SSA(aⱼ:: Vector{T})::Tuple{Int, Float64} where T<:Real
    a₀::Real = sum( aⱼ )
    j::Int = 0
    if a₀ <= 0
        println("womp womp")
        return j, 0.
    end

    r = rand(Uniform(eps(Float64),1), 2)
    τ::Float64 = log(1 / r[1]) / a₀
    rₐ₀::Float64 = r[2] * a₀

    for (i, a) ∈ enumerate( aⱼ)
        if sum( aⱼ[1:i] ) >= rₐ₀
            j = i
            break
        end
    end

    return j, τ
end



