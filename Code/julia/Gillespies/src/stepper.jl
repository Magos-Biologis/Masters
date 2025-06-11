"""
The function that actually does the iteration loop
"""
function StepIterator(
        model::String,
        steps::Int,
        x₀::Vector{Int},
        parameters::ParameterTypeSSA,
    )::StepperStruct

    varNum = length(x₀)
    parameters /= varNum # This works cause I defined how the struct divides


    ## Instantiating all the outputs
    time_array  = Vector{Real}(undef, steps)
    state_array = Array{Int}(undef, steps, varNum)
    time_array[1]    = 0
    state_array[1,:] = x₀


    ## This works because I defined the `prepropensity` to output
    ## a function with the appropriate parameters filled in
    ## I'm sure python can do it too, but it's much smoother in Julia
    ( aⱼ, v ) = instantiate_propensities(model, parameters)


    for i ∈ 2:steps
        x = state_array[i-1, :]
        (j, τ) = SSA( aⱼ(x) )

        if j == 0; final_step = i; break; end

        time_array[i] = time_array[i-1] + τ
        state_array[i, :] = x + v[j]
    end

    output = StepperStruct(time_array, state_array)
    return output
end
