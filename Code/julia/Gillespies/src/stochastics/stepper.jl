"""
The function that actually does the iteration loop
using Base: final_shred!
"""
function StepIterator(
        model::String,
        steps::Int64,
        x₀::Vector{Int64},
        parameters::ParameterType,
    )::StepperStruct

    ## Instantiating the parameters
    var_count = sum(x₀)
    parameters /= var_count # This works cause I defined how the struct divides


    ## Instantiating all the outputs
    distinct_vars = length(x₀)
    time_array  = Vector{Float64}(undef, steps)
    state_array = Array{Int32}(undef, steps, distinct_vars)
    time_array[1]    = 0
    state_array[1,:] = x₀


    ## This works because I defined the `prepropensity` to output
    ## a function with the appropriate parameters filled in
    ## I'm sure python can do it too, but it's much smoother in Julia
    ( aⱼ, v ) = instantiate_propensities(model, parameters)
    ssa_results = SSAOutput(0., 0)

    for i ∈ 2:steps
        x = state_array[i-1, :]
        ssa::SSAOutput = SSA!( ssa_results, aⱼ(x) )

        if ssa.j == 0
            final_step = i - 1
            return StepperStruct(
                            time_array[1:final_step],
                            state_array[1:final_step,:]
                       )
        end

        time_array[i] = time_array[i-1] + ssa.τ
        state_array[i, :] = x + v[ssa.j]
    end

    return StepperStruct(time_array, state_array)
end
