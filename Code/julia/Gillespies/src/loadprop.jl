function load_propensity_stuff(module_name::String)
    ## Includes the file relative to the const `PROPDIR`
    module_file::String = PROPDIR * module_name * ".jl"
    include(module_file)

end
#
# mod = getfield(Main, Symbol(module_name))
# return getfield(Symbol(module_name), :propensity)
