import PrecompileTools

PrecompileTools.@compile_workload begin
    nparm = NovelStruct(n = 100)
    dparm = DifferentialStruct(m₀ = 1)
    SimpleChemical1Par1VarAB(nparm)
    SimpleChemical1Par1VarAB(nparm; symbolic = false)
end
