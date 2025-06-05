


U::Uniform = Uniform(eps(Float64),1)




function IterationOfSSA(aⱼ:: Vector{Float64})::Tuple{Int, Float64}
    a₀::Float64 = sum( aⱼ.entries )
    if a₀ <= 0
        println("womp womp")
        return -1, 0
    end

    r::Vector{Float64} = rand(U, 2)
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
        propensities::Function,
        steps::Int,
        x₀::Vector,
        vⱼ::AbstractArray,
        parameters::SSAParameters,
    )

    m::Int = sum(x₀)
    scaled_paramters = div(parameters, m)


end
