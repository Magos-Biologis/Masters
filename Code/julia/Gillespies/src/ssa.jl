



"""
A Julia-fied version of my python implementation of the Gillespie Stochastic
Simulation Algorithm.
"""
function SSA(aⱼ:: Vector{Real})::Tuple{Int, Float64}
    a₀::Real = sum( aⱼ.entries )
    if a₀ <= 0
        println("womp womp")
        return -1, 0
    end

    r::Vector{Float64} = rand(Uniform(eps(Float64),1), 2)
    j::Int = -1
    τ::Float64 = log(1 / r[1]) / a₀
    rₐ₀::Float64 = r[2] * a₀

    for ( j, sⱼ) ∈ enumerate(cumsum( aⱼ ))

        if sⱼ>= rₐ₀
            break
        end
    end

    return j, τ
end



function StepIterater(
        model::String,
        steps::Int,
        x₀::Vector,
        parameters::ParameterTypeSSA,
        # propensities::Function,
        # vⱼ::AbstractArray,
    )

    (propensity, vⱼ)= load_propensity_stuff(model)



end

function StepSSA()
end


