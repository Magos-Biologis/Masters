"""
The function that actually does the iteration loop
"""
function StepIterator(
        model::String,
        steps::Int,
        x₀::Vector{Int},
        parameters::ParameterTypeSSA,
    )::StepperStruct

    var_count = sum(x₀)
    distinct_vars = length(x₀)
    scaled_pars = parameters / var_count # This works cause I defined how the struct divides


    ## Instantiating all the outputs
    time_array  = Vector{Float64}(undef, steps)
    state_array = Array{Int32}(undef, steps, distinct_vars)
    time_array[1]    = 0
    state_array[1,:] = x₀

    ## This works because I defined the `prepropensity` to output
    ## a function with the appropriate parameters filled in
    ## I'm sure python can do it too, but it's much smoother in Julia

    ( aⱼ, v ) = instantiate_propensities(model, scaled_pars)

    for i ∈ 2:steps
        x = state_array[i-1, :]
        (j, τ) = SSA( aⱼ(x) )

        # println(τ)
        if j == 0
            break
        end

        time_array[i] = time_array[i-1] + τ
        state_array[i, :] = x + v[j]
    end

    return StepperStruct(time_array, state_array)
end
