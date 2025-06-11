__precompile__()
module SSAExecutionModule


const BASEDIR::AbstractString = @__DIR__


using ArgParse


export ArgumentSetter,

       ParameterDict,
       ParVarStruct,

       ModelCounts,
       RateParameterPrimer

       # init




# abstract type Parameters end

struct ParameterDict{T<:Dict}
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
function ArgumentSetter()
    settings = ArgParseSettings()

    @add_arg_table! settings begin
        "model"
        help = "State which model to run"
        required = true

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

    return parse_args(settings)
end




## The model const functions
include("model_counts.jl")

## For the functions that make the parameters work
include("parameter_setting.jl")

end
