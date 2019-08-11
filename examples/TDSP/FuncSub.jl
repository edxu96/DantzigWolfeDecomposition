# Dantzig-Wolfe Reformulation and Column Generation
# Functions for Sub-Problem
# Edward J. Xu
# 2019.5.28
########################################################################################################################


mutable struct ModelSub
    mod::JuMP.Model
    mat_e
    vec_l
    mat_d
    vec_q
    vec_x
end


function setVariableSub(modSub::JuMP.Model, mat_d, vec_q, vecIndexBin)
    (m_d, n_d) = size(mat_d)
    @variable(modSub, vec_x[1:n_d] >= 0, Int)  # , Bin
    # @constraint(modSub, vec_x[vecIndexBin] .<= 1)
    # objective is not here. We define once dual variables become known
    @objective(modSub, Max, 0)
    # remember to change "<=" if your sub-problem uses a different type of constraints!
    @constraint(modSub, cons[i = 1: m_d], sum(mat_d[i, j] * vec_x[j] for j = 1:n_d) <= vec_q[i])
    return vec_x
end


function setBranch(vecModelSub, whiSub, whiVar, whiBranch)
    numXInSub = size(vecModelSub[1].mat_d)[2]
    vecBranch = zeros(1, numXInSub)
    if whiBranch == 1
        vecBranch[whiVar] = -1
        vecModelSub[whiSub].vec_q = vcat(vecModelSub[whiSub].vec_q, [-1])
    elseif whiBranch == 0
        vecBranch[whiVar] = 1
        vecModelSub[whiSub].vec_q = vcat(vecModelSub[whiSub].vec_q, [0])
    end
    vecModelSub[whiSub].mat_d = vcat(vecModelSub[whiSub].mat_d, vecBranch)
    println("$(whiVar)-th variable in $(whiSub)-th sub-problem is branched to $(whiBranch).")
    return vecModelSub
end


function setModelSub(
    numSub, vecVecIndexSub, m_mat_h, vec_c, mat_a, vec_b, vecRowMatD, gurobi_env, vecVecIndexBinInSub,
    wheBranch, vecWhiSub, vecWhiVar, vecWhiBranch
    )
    vecModelSub = Vector{ModelSub}(undef, numSub)
    for k = 1: numSub
        vecIndexCol = vecVecIndexSub[k]
        if k == 1
            vecIndexRow = collect(convert(Int64, (m_mat_h + 1)):convert(Int64, (m_mat_h + vecRowMatD[1])))
        else
            vecIndexRow = collect(convert(Int64, (m_mat_h + sum(vecRowMatD[1:(k-1)]) + 1)):
                convert(Int64, (m_mat_h + sum(vecRowMatD[1:k]))))
        end
        # println("typeof(vecIndexCol) = $(typeof(vecIndexCol)), typeof(vecIndexRow) = $(typeof(vecIndexRow)).")
        vecModelSub[k] = ModelSub(
            Model(solver = GurobiSolver(OutputFlag = 0, gurobi_env)),                          # mod
            mat_a[1: m_mat_h, vecIndexCol],                                              # mat_e,
            vec_c[vecIndexCol],                                                          # vec_l
            hcat(mat_a[vecIndexRow, vecIndexCol]),   # mat_d
            vec_b[vecIndexRow, 1],                         # vec_q
            0                                                                                  # vec_x
            )
    end
    # Set Branch
    if wheBranch
        for i = 1: length(vecWhiVar)
            vecModelSub = setBranch(vecModelSub, vecWhiSub[i], vecWhiVar[i], vecWhiBranch[i])
        end
    end
    for k = 1: numSub
        vecModelSub[k].vec_x = setVariableSub(vecModelSub[k].mod, vecModelSub[k].mat_d, vecModelSub[k].vec_q,
            vecVecIndexBinInSub[k])
    end
    return vecModelSub
end


function solveSub(modSub::ModelSub, vec_pi, kappa)
    (m, n) = size(modSub.mat_e)
    piA0 = vec_pi * modSub.mat_e
    @objective(modSub.mod, Max, sum(modSub.vec_l[j] * modSub.vec_x[j] for j = 1: n) -
        sum(piA0[1, j] * modSub.vec_x[j] for j = 1: n) - kappa)
    status = solve(modSub.mod)
    if status != :Optimal
        throw("Error: Non-optimal sub-problem status")
    end
    costReduce = getobjectivevalue(modSub.mod)
    vec_x_result = getvalue(modSub.vec_x)
    return (costReduce, vec_x_result)
end
