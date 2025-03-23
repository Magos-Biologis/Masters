#using Random: UniformT
# const output = IOBuffer()
# using REPL
# const out_terminal = REPL.Terminals.TerminalBuffer(output)
# const basic_repl = REPL.BasicREPL(out_terminal)
# const basic_display = REPL.REPLDisplay(basic_repl)

# Base.pushdisplay(basic_display)

### Packages

using Plots;
gr();

using StatsPlots
using Distributions

using LaTeXStrings
using Random

#using DifferentialEquations


### variables

file_path = "figs/test.png"

## base
w::Vector{Float64} = [0; 0]
k::Vector{Float64} = [0; 0]
n::Vector{Float64} = [0; 0]

q::Vector{Float64} = [0; 0]

m_0::Float64 = 0

## subs
l::Vector{Float64} = [0; 0]
ω::Vector{Float64} = [0; 0]


t_end = 5000

## vals
k[1] = 2
k[2] = 2

# Population cap (purely aesthetic if n₁ = n₂)
n[1] = 100
n[2] = 100

# w[1] = 0.015
# w[2] = 0.035

w[1] = 0.7
w[2] = 0.4

q[1] = 0.999
q[2] = 0.8


l[1] = k[1] / n[1]
l[2] = k[2] / n[1]

ω[1] = k[1] - w[1]
ω[2] = k[2] - w[2]

## initial values
c₁₀ = 0.99
c₂₀ = 1 - c₁₀

c::Vector{Float64} = [c₁₀; c₂₀]

state_vector = c

### Chain

test_matrix::Matrix{Float64} = [1-w[1] w[1]; w[2] 1-w[2]]
#test_matrix::Matrix{Float64} = [(ω[1]-l[1]) w[1]; w[2] (ω[2]-l[2])]

# x_array = Vector{Union{Int,Missing}}(undef, t_end)
# y_array = Vector{Union{Int,Missing}}(undef, t_end)

x_array::Vector{Int} = zeros(t_end)
y_array::Vector{Int} = zeros(t_end)
t_array::Vector{Int} = 1:1:t_end

x_array[1] = 99
y_array[1] = 1

# println(t_array)
# exit()

# append!(x_array, 0)
# append!(y_array, 0)
# append!(t_array, 0)

is_c1 = 1
is_c2 = 2


mutable struct StateTracking{Int}
    state::Int
end
markov = StateTracking(is_c1)



# time_modifier = [1]

Random.seed!(190622)

U = Uniform(0, 1)
u = rand(U, 2, length(t_array[1:end-1]))
#
# print(u)
# exit()
# println(t_array)
# exit()

#while t_array[end] <= t_end
for tau in t_array[1:end-1]
    local upsilon = tau + 1

    x_array[upsilon] += x_array[tau]
    y_array[upsilon] += y_array[tau]

    if markov.state == is_c1
        if u[1, tau] < test_matrix[1, 1]
            x_array[upsilon] += 1
            # y_array[upsilon] = y_array[tau]
        else
            # x_array[upsilon] = x_array[tau]
            y_array[upsilon] += 1

            markov.state = is_c2
        end
    end

    if markov.state == is_c2
        if u[2, tau] < test_matrix[2, 2]
            # x_array[upsilon] = x_array[tau]
            y_array[upsilon] += 1
        else
            x_array[upsilon] += 1
            # y_array[upsilon] = y_array[tau]

            markov.state = is_c1
        end
    end

    ## End of chain, increments at end of loop so arrays work properly
    #tau[1] += 1
end


sums = x_array + y_array
x::Vector{Float64} = x_array ./ sums
y::Vector{Float64} = y_array ./ sums

# println(sums)
# exit()



# print(x_array)
# exit()





Plots.default(linewidth=2, ylims=(0, 1))

test_plot = plot(t_array, x; label=(L"x"))
plot!(t_array, y; label=(L"y"))
#plot!(t_array, y)
# ylims!(0, 1)
xlims!(0, Inf)

savefig(test_plot, file_path)

println("done")
exit()
