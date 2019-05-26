

function doDantzigWolfeDecomp(vecModelSub, modMas, vecConsRef, vecConsConvex,
    vecLambda, vecStrNameVar, epsilon, whePrint, indexSub)
    numSub = length(vecModelSub)
    extremePoints = [[]]
    extremePointForSub = [-1]
    done = false
    iter = 1
    while !done
        (vec_pi, vec_kappa, obj_master) = solveMas(modMas, vecConsRef, vecConsConvex)
        done = true
        ## Print the result
        if whePrint
            println("vec_pi = $(vec_pi)\n",
                    "vec_kappa = $(vec_kappa)")
        end
        ## Column Generation
        costReduceBest = -1
        for k = 1: numSub
            (costReduce, vec_x_result) = solveSub(vecModelSub[k], vec_pi, vec_kappa[k])
            if whePrint
                println("Reduced cost of $(k)-th sub model is $costReduce.")
            end
            if costReduce > costReduceBest
                costReduceBest = costReduce
            end
            # remember to look for negative reduced cost if we are minimizing.
            if costReduce > epsilon
                lambdaNew = addColToMaster(modMas, vecModelSub[k], vec_x_result, vecConsRef, vecConsConvex[k])
                push!(vecLambda, lambdaNew)
                push!(extremePoints, vec_x_result)
                push!(extremePointForSub, k)
                done = false
            end
        end
        ##
        iter += 1
        println("best reduced cost is $costReduceBest")
        println("$(iter)-th iteration. obj_master = $(obj_master).\n",
                "--------------------------------------------------------------------------------")
    end
    ## 5,  Print the Result
    println("################################################################################\n", ######################
            "Optimization Done After $(iter-1) Iterations.\n",
            "objMaster = $(getobjectivevalue(modMas))\n",
            "################################################################################")
    # compute values of original variables
    origVarValSub = []
    for s = 1: length(indexSub)
        push!(origVarValSub, zeros(length(indexSub[s])))
    end
    vec_lambda_result = getvalue(vecLambda)
    for p = 1: length(vec_lambda_result)
        if vec_lambda_result[p] > 0.000001
            println("lambda_$p = ", vec_lambda_result[p], ", sub = $(extremePointForSub[p]), \n",
                    "extreme point = $(extremePoints[p]).")
            origVarValSub[extremePointForSub[p]] += vec_lambda_result[p] * extremePoints[p]
        end
    end
    for s = 1: length(indexSub)
        #println("var val for sub problem $s: $(origVarValSub[s])")
        for t = 1: length(origVarValSub[s])
            if abs(origVarValSub[s][t]) > 0.000001
                println("$(vecStrNameVar[indexSub[s][t]]) = $(origVarValSub[s][t])")
            end
        end
    end
end
