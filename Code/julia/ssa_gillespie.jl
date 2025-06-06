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

### The valid arguments
model_path = getfield(Gillespies, PROPDIR)
files = readdir(model_path)

file_set = Set{String}()
for file in files
    file = replace(file, ".jl" => "")
    push!(file_set, file)
end

function parse_parameters(param_str::String)
    matches = eachmatch(r"(\w+[^=]*)=\s?([^ ,]*)", param_str)
    return Dict(
        Symbol(replace(m.captures[1], ["-" => "_", "'" => "p", " " => ""])) => parse(Float64, m.captures[2])
        for m in matches
    )
end

function parse_args()
    s = ArgParseSettings()

    @add_arg_table s begin
        "model"
        help = "State which model to run"
        arg_type = String
        choices = collect(file_set)

        "--repeats", "-r"
        help = "Number of repetitions"
        arg_type = Int
        default = 1

        "--no-save", "-ns"
        help = "Flag to disable saving"
        action = :store_false
        dest_name = :save

        "--test-parameters"
        help = "Prints the parameter dictionary and exits"
        action = :store_true
        dest_name = :test

        "--steps", "-s"
        help = "Number of steps"
        arg_type = Int
        default = 10_000

        "--size", "-n"
        help = "Total units in system"
        arg_type = Int

        "--initial_conds", "-ic"
        help = "Initial conditions"
        nargs = '+'
        arg_type = Int
        default = [99, 0, 0]

        "--parameters", "-p"
        help = "Parameter string like 'k=1.0, alpha=0.1'"
        arg_type = parse_parameters
        default = Dict(:k => 1.0)  # replace with real default if needed
    end

    return parse_args(s)
end

args = parse_args()

if get(args, :test, false)
    println(args[:parameters])
    exit()
end
