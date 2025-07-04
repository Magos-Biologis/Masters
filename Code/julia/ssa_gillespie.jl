# ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0    # I don't think this is needed if
# import Pkg                              # I do `julia --project` when
# Pkg.activate(@__DIR__)                  # running the file
# ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 1    #


# pushfirst!(LOAD_PATH, joinpath(BASEDIR, "Gillespies/"))
# pushfirst!(LOAD_PATH, joinpath(BASEDIR, "SSAExecutionModule/"))

const BASEDIR::AbstractString = @__DIR__
const DATASTORE::AbstractString = get(ENV, "THESIS_DATA_PATH", "")

using ArgParse
using JSON
using NPZ

using Random

using Gillespies #-> Custom local module
using SSAExecutionModule #-> Custom local module


parsed_args = ArgumentSetter()



model = parsed_args["model"]
is_ode = occursin(r"((ODE)|(ode))", model)

params = parsed_args["parameters"]

steps = parsed_args["steps"]
initial_condition = parsed_args["init"]

repeats = parsed_args["repeats"]

save_file = parsed_args["save"]


file_name = "ssa"
file_name *= "M" * model
file_name *= "L" * "julia"


# TODO: Now that the parameters are parsed properly, I want to make sure they
# get stored appropriatly, like k₄ → k⁺[4], and k₋₂ → into k⁻[2].
# Additionally that all the parameters go into the eventual 'metadata' .json
# for .npz

const relevant_counts::ParVarStruct = ModelCounts(model)
initial::Vector{<:Integer} = initial_condition[1:relevant_counts.var]

particle_count = sum(initial)
rate_vectors = RateParameterPrimer(params, is_ode, relevant_counts)
if is_ode
    parameters = DifferentialStruct(;
        m₀=get(params.all, "m0", 0),
        n = rate_vectors[1],
        k = rate_vectors[2],
        w = rate_vectors[3],
        q = rate_vectors[4],
    )
else
    parameters = NovelStruct(;
        n = get(params.all, "n", 1),
        b = get(params.all, "b", 0),
        k⁺= rate_vectors[1],
        k⁻= rate_vectors[2],
    )
end


"""
Setting the metadata json
"""
metadata_dict = Dict{String, Any}()
merge!(metadata_dict,
       Dict(
            "data_source" => "stochastic simulation algorithm",
            "model_name" => model,
            "number_of_variables" => relevant_counts.var,
            "initial_condition" => initial,
            "number_of_particles" => sum(initial),
            "parameters" => params.all
           )
      )


for i::Int64 in 1:repeats
    save_name = deepcopy(file_name)
    metadata_cycle = deepcopy(metadata_dict)

    t1 = Float64(time())
    results = StepIterator(model, steps, initial, parameters)
    t2 = Float64(time())

    Δt::Float64 = t2 - t1
    steps_taken::Int64 = length(results.time)

    # println(Δt)


    epoch = t1
    save_name *= "T" * replace("$epoch", "." => "", "e9" => "")

    update_dict = Dict(
                       "steps_taken" => steps_taken,
                       "date" => epoch,
                       "runtime" => Δt,
                       "run" => i,
                      )

    merge!(metadata_cycle, update_dict)
    metadata_json = json(metadata_cycle)
    metadata_bytes = collect(codeunits(metadata_json))

    to_write = Dict(
                    "time"     => results.time,
                    "states"   => results.states,
                    "metadata" => metadata_bytes
                   )


    full_file_path = joinpath(DATASTORE, save_name)
    if save_file
        npzwrite(
             full_file_path * ".npz",
             to_write
            )
        println("saved as '$save_name'")
    else
        println("file name would be is '"*"$save_name"*".npz"*"'")
    end

end

exit()
