using ArgParse
using JSON
using Dates
using Printf
using Random
using Glob
using Pkg

# Custom deepcopy import
import Base.deepcopy as dc

using Gillespies

# Required for environment variables
figure_env = get(ENV, "FPE_FIGURE_ENV", "")
data_env = get(ENV, "THESIS_DATA_PATH", "")
code_env = get(ENV, "THESIS_CODE", "")

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
    positive = filter(x -> x.captures[2] == nothing, matches)
    negative = filter(x -> x.captures[2] != nothing, matches)
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
    arg_type = Int
    default = 1

    "--steps", "-s"
    help = "Number of steps"
    arg_type = Int
    default = 10_000

    "--size", "-n"
    help = "Total units in system"
    arg_type = Int

    "--initial_conds", "-i"
    help = "Initial conditions"
    nargs = '+'
    default = [99, 0, 0]

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


# TODO: Now that the parameters are parsed properly, I want to make sure they
# get stored appropriatly, like k₄ → k⁺[4], and k₋₂ → into k⁻[2].
# Additionally that all the parameters go into the eventual 'metadata' .json
# for .npz



exit()




if get(args, :test, false)
    println(args[:parameters])
    exit()
end
