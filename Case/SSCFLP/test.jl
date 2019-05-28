# Dantzig-Wolfe Decomposition
# Edward J. Xu
# 2019.5.28
# include("$(homedir())/Desktop/9-SSCFLP/test.jl")
########################################################################################################################
push!(LOAD_PATH, "$(homedir())/Desktop/9-SSCFLP")
cd("$(homedir())/Desktop/9-SSCFLP")
using JuMP
# using CPLEX
using Gurobi
using LinearAlgebra
# using DantzigWolfeOptim
include("FuncSub.jl")
include("FuncMas.jl")
include("FuncOptim.jl")
include("FuncData.jl")
include("DantzigWolfeOptim.jl")
########################################################################################################################


function main()
    ## 1,  Data
    (vec_c, mat_a, vec_b, indexSub, numSub, numXInSub, m_mat_h, vecRowMatD) = getData(1)
    ## 2,  Optim
    wheBranch = true
    vecWhiSub = [1, 1]
    vecWhiVar = [12, 6]
    vecWhiBranch = [1, 1]
    doDantzigWolfeOptim(vec_c, mat_a, vec_b, indexSub, numSub, numXInSub,
        m_mat_h, vecRowMatD, wheBranch, vecWhiSub, vecWhiVar, vecWhiBranch
        )
end


main()
