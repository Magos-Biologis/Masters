"""
This function is used to load the individual modules for the
simulation and extract the relevant 'fields'.
"""
function load_propensity_stuff(module_name::AbstractString)::Tuple{Function, AbstractArray}
    module_file_path::AbstractString = joinpath(PROPDIR, module_name * ".jl")

    include(module_file_path)
    symbol_name::Symbol = Symbol(module_name)
    mod = getfield(Gillespies, symbol_name)

    propensity = getfield(mod, :propensity)
    transitions = getfield(mod, :transitions)

    return propensity, transitions
end


"""
Includes the file relative to the const `PROPDIR`

I originally tried to make this load a module in the function,
similar to how python can do an `import_module` → `getattribute`.
However, Julia does not support modules at anything less that the
top level, so I had to simplify it to just include the file.

I take it back, Julia does support something similar, it was just that
I kept using a number for the start of the module name, which
is unsupported.
"""
# println(propensity, transitions)

