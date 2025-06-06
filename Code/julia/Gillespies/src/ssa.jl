"""
A Julia-fied version of my python implementation of the Gillespie Stochastic
Simulation Algorithm.
"""
function SSA(aⱼ:: Vector{T})::Tuple{Int, Float64} where T<:Real
    a₀::Real = sum( aⱼ )
    if a₀ <= 0
        println("womp womp")
        return 0, 0
    end

    r::Vector{Float64} = rand(Uniform(eps(Float64),1), 2)
    j::Int = 0
    τ::Float64 = log(1 / r[1]) / a₀
    rₐ₀::Float64 = r[2] * a₀

    for ( i, sⱼ) ∈ enumerate(cumsum( aⱼ ))
        if sⱼ>= rₐ₀
            j += i
            break
        end
    end

    return j, τ
end



