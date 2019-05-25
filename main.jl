# Dantzig-Wolfe Reformulation and Column Generation
# Edward J. Xu
# 2019.5.1
########################################################################################################################
# push!(LOAD_PATH, "$(homedir())/Documents/Github/DantzigWolfeOptim")
# cd("$(homedir())/Documents/Github/DantzigWolfeOptim")
using JuMP
# using CPLEX
using Gurobi
using LinearAlgebra
include("FuncSub.jl")
include("FuncMas.jl")
########################################################################################################################


function addColToMaster(mod_mas::JuMP.Model, modSub::ModelSub, vec_x_result, vec_consRef, consConvex)
    (m, n) = size(modSub.mat_e)
    A0x = modSub.mat_e * vec_x_result
    vec_consTouched = ConstraintRef[]
    vals = Float64[]
    for i = 1: m
        # only insert non-zero elements (this saves memory and may make the master problem easier to solve)
        if A0x[i] != 0
            push!(vec_consTouched, vec_consRef[i])
            push!(vals, A0x[i])
        end
    end
    # add variable to convexity constraint.
    push!(vec_consTouched, consConvex)
    push!(vals, 1)
    objCoef = sum(modSub.vec_l[j] * vec_x_result[j] for j = 1: n)
    println("objCoef = $objCoef.")
    @variable(
        mod_mas,
        lambdaNew >= 0,                     # New variable to be added
        objective = objCoef,                # cost coefficient of new varaible in the objective
        inconstraints = vec_consTouched,    # constraints to be modified
        coefficients = vals                 # the coefficients of the variable in those constraints
        )
    return lambdaNew
end


function doDantzigWolfeReform(vec_c, mat_a, vec_b, index_sub, num_sub, num_x_sub, num_border, vecStr_nameVar, epsilon)
    ## 2,  Set up all the sub models
    vecModelSub = setModelSub(num_sub, num_border, vec_c, mat_a, vec_b)
    ## 3,  Set Up Master Models
    (mod_mas, vec_consRef, vec_consConvex, vec_lambda) = setModelMas(vecModelSub[1].mat_e, vec_b, num_sub, num_border)
    ## 4,  Begin Iteration
    println("#### Begin Iteration ###########################################################") ########################
    extremePoints = [[]]
    extremePointForSub = [-1]
    done = false
    iter = 1
    while !done
        (vec_pi, vec_kappa, obj_master) = solveMas(mod_mas, vec_consRef, vec_consConvex)
        done = true
        ## Print the result
        println("vec_pi = $(vec_pi)\n",
                "vec_kappa = $(vec_kappa)")
        println("--------------------------------------------------------------------------------\n",
                "$(iter)-th iteration. obj_master = $(obj_master).\n",
                "--------------------------------------------------------------------------------")
        ## Column Generation
        costReduceBest = -1
        for k = 1: num_sub
            (costReduce, vec_x_result) = solveSub(vecModelSub[k], vec_pi, vec_kappa[k])
            println("Reduced cost of $(k)-th sub model is $costReduce.")
            if costReduce > costReduceBest
                costReduceBest = costReduce
            end
            # remember to look for negative reduced cost if we are minimizing.
            if costReduce > epsilon
                lambdaNew = addColToMaster(mod_mas, vecModelSub[k], vec_x_result, vec_consRef, vec_consConvex[k])
                push!(vec_lambda, lambdaNew)
                push!(extremePoints, vec_x_result)
                push!(extremePointForSub, k)
                done = false
            end
        end
        ##
        iter += 1
        println("best reduced cost is $costReduceBest")
    end
    ## 5,  Print the Result
    println("################################################################################\n", ######################
            "Optimization Done After $(iter-1) Iterations.\n",
            "obj_master = $(getobjectivevalue(mod_mas))\n",
            "################################################################################")
    # compute values of original variables
    origVarValSub = []
    for s = 1: length(index_sub)
        push!(origVarValSub, zeros(length(index_sub[s])))
    end
    vec_lambda_result = getvalue(vec_lambda)
    for p = 1: length(vec_lambda_result)
        if vec_lambda_result[p] > 0.000001
            println("lambda_$p = ", vec_lambda_result[p], ", sub = $(extremePointForSub[p]), \n",
                    "extreme point = $(extremePoints[p]).")
            origVarValSub[extremePointForSub[p]] += vec_lambda_result[p] * extremePoints[p]
        end
    end
    for s = 1: length(index_sub)
        #println("var val for sub problem $s: $(origVarValSub[s])")
        for t = 1: length(origVarValSub[s])
            if abs(origVarValSub[s][t]) > 0.000001
                println("$(vecStr_nameVar[index_sub[s][t]]) = $(origVarValSub[s][t])")
            end
        end
    end
end


function main()
    timeStart = time()
    epsilon = 0.00001
    global gurobi_env = Gurobi.Env()
    ## 1,  Data
    include("Data_GAP_2.jl")
    ## 2,  Optimization
    doDantzigWolfeReform(vec_c, mat_a, vec_b, index_sub, num_sub, num_x_sub, num_border, vecStr_nameVar, epsilon)
    println("################################################################################\n",
            "Elapsed time is $(time() - timeStart) seconds.\n",
            "################################################################################")
end


main()
