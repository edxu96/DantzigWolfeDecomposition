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


function main(strFileName::String, series::Int64, whePrint::Bool)
    strFileName = "data/gapa.txt"
    series = 1
    timeStart = time()
    epsilon = 0.00001
    global gurobi_env = Gurobi.Env()
    println("#### 1/5,  Read Data ###########################################################") ########################
    (w, p, c) = readData(strFileName, series)
    (modAll, numXPerSub, blocks) = setGAPModel(w, p, c)
    println("#### 2/5,  Get Index ###########################################################") ########################
    (mat_a, vec_b, vec_c, vecSense) = getMat(modAll)
    (indexMas, indexSub) = getIndex(mat_a, blocks)
    println("#### 3/5,  Set vecModelSub #####################################################") ########################
    vecModelSub = setModelSub(mat_a, vec_b, vec_c, vecSense, indexMas, blocks, indexSub)
    numSub = length(vecModelSub)
    numQ = size(vecModelSub[1].mat_e)[1]  # [num of constraints in matA_0]
    # vecModelSub = setBranch(vecModelSub, 1, 2)
    for k = 1:numSub
        vecModelSub[k].vec_x = setVariableSub(vecModelSub[k].mod, vecModelSub[k].mat_d, vecModelSub[k].vec_q,
            vecModelSub[k].vecSense)
    end
    vecStrNameVar = getStrNameVar(numSub, numXPerSub)  # define variable names
    println("#### 4/5,  Set modelMas ########################################################") ########################
    ## Parameters for stabilization
    vecDualGuessPi = 100 * ones(Float64, numQ)
    vecDualGuessKappa = 0 * ones(Float64, numSub)
    dualPen = 10
    dualPenMult = 0.1
    dualPenThreshold = 0.01 - 1e-5
    ##
    vecSenseP = deepcopy(vecSense[collect(indexMas)])
    vecP = deepcopy(vec_b[collect(indexMas)])
    (modMas, vecConsRef, vecConsConvex, vecLambda, vecMuMinus, vecMuPlus, vecMuMinusConv, vecMuPlusConv) =
        setModelMas(numQ, vecP, numSub, vecSenseP, dualPen, vecDualGuessPi, vecDualGuessKappa)
    vecObjCoef = [-1000.0]
    ## 5,  Optimization
    println("#### 5/5,  Begin Optim #########################################################") ########################
    doDantzigWolfeDecomp(vecModelSub, modMas, vecConsRef, vecConsConvex, vecLambda, vecStrNameVar, epsilon, whePrint,
        indexSub, dualPen, dualPenMult, dualPenThreshold, vecDualGuessPi, vecDualGuessKappa, vecMuMinus, vecMuPlus,
        vecMuMinusConv, vecMuPlusConv, vecObjCoef)
    println("################################################################################\n",
            "Elapsed time is $(time() - timeStart) seconds.\n",
            "################################################################################")
end


@enum Sense leq = 1 geq = 2 eq = 3
main("data/gapa.txt", 1, false)
