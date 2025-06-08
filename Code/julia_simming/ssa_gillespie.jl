

using ArgParse
using JSON
using NPZ


using Gillespies


const BASEDIR = @__DIR__

# Required for environment variables
figure_env = get(ENV, "FPE_FIGURE_ENV", "")
data_env = get(ENV, "THESIS_DATA_PATH", "")
code_env = get(ENV, "THESIS_CODE", "")

# println(data_env)
# exit()


# directories = ["pypacks", "gillespie", "propensities", "aj"]
# ### The valid arguments
# my_module = getfield(Main, :Gillespies)
# model_path = getfield(my_module, :PROPDIR)
# files = readdir(model_path)

# file_set = Set{String}()
# for file in files
#     file = replace(file, ".jl" => "")
#     push!(file_set, file)
# end
#

abstract type Parameters end

struct ParameterDict{T<:Dict} <: Parameters
    positive::T
    negative::T
    all::T
end



function directional_rate(matches::Vector{RegexMatch}) :: Tuple{Vector{RegexMatch}, Vector{RegexMatch}}
    positive = filter(x ->  isnothing(x.captures[2]), matches)
    negative = filter(x -> ~isnothing(x.captures[2]), matches)
    return positive, negative
end

function parse_parameters(parameter_string::AbstractString) :: ParameterDict
    par_regex = r"(\w+(-)?[^=]*)=\s?([^ ,]*)"
    replaces = ["-" => "_", "'" => "p", " " => ""]
    pairing(match::RegexMatch) = replace(match.captures[1], replaces...) => parse(Float64, match.captures[3])

    matches = collect(eachmatch(par_regex, parameter_string))

    pos, neg = directional_rate(matches)
    pos_pairs, neg_pairs, pairs = map(pairing, pos), map(pairing, neg), map(pairing, matches)
    return ParameterDict(Dict(pos_pairs), Dict(neg_pairs), Dict(pairs))
end

## Overloading the argparse typing
function ArgParse.parse_item(::Type{ParameterDict}, x::AbstractString)
    return parse_parameters(x)
end



"""
    Parsing the arguments given from the command line, trying to be as similar
    as possible to the python ones.
"""
settings = ArgParseSettings()

@add_arg_table! settings begin
    "model"
    help = "State which model to run"
    required = true
    # action = :command

    "--repeats", "-r"
    help = "Number of repetitions"
    arg_type = Int64
    default = 1

    "--steps", "-s"
    help = "Number of steps"
    arg_type = Int64
    default = 10_000

    "--size", "-n"
    help = "Total units in system"
    arg_type = Int

    "--initial-conds", "-i"
    help = "Initial conditions"
    nargs = '+'
    default = [99; 0; 0]
    dest_name = "init"
    arg_type = Int64

    "--parameters", "-p"
    help = "Parameter string like 'k=1.0, alpha=0.1'"
    arg_type = ParameterDict
    # default = Dict(:k => 1.0)  # replace with real default if needed
    # arg_type = parse_parameters

    ## Boolean flags
    "--no-save"
    help = "Flag to disable saving"
    action = :store_false
    dest_name = "save"

    "--test-parameters"
    help = "Prints the parameter dictionary and exits"
    action = :store_true
    dest_name = "test"
end

parsed_args = parse_args(settings)

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


## The model const functions
include(joinpath(BASEDIR, "model_counts.jl"))


"""
This function modifies a given vector in place using the Dict.
By taking the value of the 'subscript' of a parameter, it places that one
in the (hopefully) correct place the vector
"""
function parameter_assignment!(par_vec::Vector, parameters::Dict)#, is_ode::Bool)
    for (k, v) ∈ parameters
        entry = match(r"(\d+)", String(k)).captures[1]
        par_vec[ parse(Int, entry) ] = v
    end
end


"""
This function constructs the missing entries in the dict
"""
function condic(dict::Dict, param::AbstractString, len::Int)
    filtered = Dict()
    for i ∈ 1:len
        pos_key = Symbol(param * "$i")
        neg_key = Symbol(param * "_$i")
        filtered[pos_key] = get(dict, pos_key, 1.)
        filtered[neg_key] = get(dict, neg_key, 1.)
    end
    return filtered
end


"""
This function filters out the relevant `Dict` entries
"""
function condic(dict::Dict, param::AbstractString, len::Int)
    filtered = Dict()
    for i ∈ 1:len
        pos_key = Symbol(param * "$i")
        neg_key = Symbol(param * "_$i")
        filtered[pos_key] = get(dict, pos_key, 1.)
        filtered[neg_key] = get(dict, neg_key, 1.)
    end
    return filtered
end


"""
This function makes sure that the parameters fed into the simulators are
consistently sized and sorted.
It could probably be made more elegant with some for-loops
"""
function rate_parameter_priming(P::ParameterDict, ode::Bool)
    if ode
        par = P.all

        n_vec = Vector{Float64}(undef, 2); n_vec .= 1
        k_vec = Vector{Float64}(undef, 2); k_vec .= 1
        w_vec = Vector{Float64}(undef, 2); w_vec .= 1
        q_vec = Vector{Float64}(undef, 2); q_vec .= 1

        parameter_assignment!( n_vec, filtdic(par, "n", 2) )
        parameter_assignment!( k_vec, filtdic(par, "k", 2) )
        parameter_assignment!( w_vec, filtdic(par, "w", 2) )
        parameter_assignment!( q_vec, filtdic(par, "q", 2) )

        return n_vec, k_vec, w_vec, q_vec
    else
        pos = P.positive
        neg = P.negative

        pos_vec = Vector{Float64}(undef, relevant_counts.par); pos_vec .= 1
        neg_vec = Vector{Float64}(undef, relevant_counts.par); neg_vec .= 1

        parameter_assignment!(pos_vec, pos)
        parameter_assignment!(neg_vec, neg)

        return pos_vec, neg_vec
    end
end


rate_vectors = rate_parameter_priming(params, is_ode)
if is_ode
    parameters = DifferentialStructSSA(;
        m₀=get(params.all, "m0", 0),
        n = rate_vectors[1],
        k = rate_vectors[2],
        w = rate_vectors[3],
        q = rate_vectors[4],
    )
else
    parameters = NovelStructSSA(;
        n = get(params.all, "n", 10),
        b = get(params.all, "b",  0),
        k⁺= rate_vectors[1],
        k⁻= rate_vectors[2],
    )
end



"""
Setting the metadata json
"""

init = initial_condition[1:relevant_counts.var]

metadata_dict = Dict{String, Any}()
merge!(metadata_dict,
       Dict(
            "data_source" => "stochastic simulation algorithm",
            "model_name" => model,
            "number_of_variables" => relevant_counts.var,
            "initial_condition" => init,
            "number_of_particles" => sum(init),
            "parameters" => params.all
           )
      )

# println(metadata_dict)
# exit()




for i::Int64 in 1:repeats
    save_name = deepcopy(file_name)
    metadata_cycle = deepcopy(metadata_dict)

    t1 = Float64(time())
    results = StepIterator(model, steps, init, parameters)
    t2 = Float64(time())

    # println(typeof(results))
    # exit()

    steps_taken = length(results.time)

    epoch = t1
    save_name *= "T" * replace("$epoch", "." => "", "e9" => "")

    update_dict = Dict(
                       "steps_taken" => steps_taken,
                       "date" => epoch,
                       "runtime" => t2 - t1,
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

    # println(to_write["time"])
    # println(to_write["states"])
    exit()

    full_file_path = joinpath(data_env, save_name)
    if save_file
        npzwrite(
             full_file_path * ".npz",
             to_write
            )
        println("saved as '$save_name'")
    else
        println("file name is '$save_name'")
    end

end

exit()
