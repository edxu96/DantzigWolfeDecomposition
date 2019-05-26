# Dantzig-Wolfe Reformulation and Column Generation
# Edward J. Xu
# 2019.5.26
########################################################################################################################
# include("$(homedir())/Documents/Github/DantzigWolfeDecomposition/main.jl")
push!(LOAD_PATH, "$(homedir())/Documents/Github/DantzigWolfeDecomposition")
cd("$(homedir())/Documents/Github/DantzigWolfeDecomposition")
using JuMP
# using CPLEX
using Gurobi
using LinearAlgebra
using MathProgBase
using Random
using SparseArrays
include("FuncSub.jl")
include("FuncMas.jl")
include("FuncDW.jl")
include("FuncData.jl")
include("FuncProcess.jl")
include("FuncStab.jl")
########################################################################################################################


function main(strFileName::String, series::Int64)
    timeStart = time()
    global gurobi_env = Gurobi.Env()
    println("#### 1/5,  Prepare Data ########################################################") ########################
    (mat_a, vec_b, vec_c, vecSenseAll, indexMas, blocks, indexSub, numXPerSub) = setGeneralAssignProbDate(
        strFileName, series)
    ## Set parameter and start DW-Decomposition
    dualPen = 10
    dualPenMult = 0.1
    dualPenThreshold = 0.01 - 1e-5
    epsilon = 0.00001
    whePrint = false
    doDWDecomp(
        mat_a, vec_b, vec_c,                                  # Data in LP Problem
        vecSenseAll, indexMas, blocks, indexSub, numXPerSub,  # Data for DW-Decomp
        dualPen, dualPenMult, dualPenThreshold,               # Para for Stable
        epsilon, whePrint                                     # Control Para
        )
    println("Elapsed time is $(time() - timeStart) seconds.\n",
            "################################################################################")
end


@enum Sense leq = 1 geq = 2 eq = 3
main("Data/gapa.txt", 1)
