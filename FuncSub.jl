# Dantzig-Wolfe Reformulation and Column Generation
# Functions for Sub-Problem
# Edward J. Xu
# 2019.5.25
########################################################################################################################


mutable struct ModelSub
    mod::JuMP.Model
    mat_e
    vec_l
    mat_d
    vec_q
    vec_x
    vecSense
end


function setVariableSub(modSub::JuMP.Model, mat_d, vec_q, vecSense)
    (m_d, n_d) = size(mat_d)
    @variable(modSub, vec_x[1: n_d] >= 0, Bin)
    @objective(modSub, Max, 0)
    # @constraint(modSub, cons[i = 1: m_d], sum(mat_d[i, j] * vec_x[j] for j = 1: n_d) <= vec_q[i])
    for i = 1:m_d
        if vecSense[i] == leq
            @constraint(modSub, sum(mat_d[i, j] * vec_x[j] for j = 1:n_d) <= vec_q[i])
        elseif vecSense[i] == geq
            @constraint(modSub, sum(mat_d[i, j] * vec_x[j] for j = 1:n_d) >= vec_q[i])
        else
            @constraint(modSub, sum(mat_d[i, j] * vec_x[j] for j = 1:n_d) == vec_q[i])
        end
    end
    return vec_x
end


function setBranch(vecModelSub, i, j)
    vec_new = zeros(1, 15)
    vec_new[j] = -1
    vecModelSub[i].mat_d = vcat(vecModelSub[i].mat_d, vec_new)
    vecModelSub[i].vec_q = vcat(vecModelSub[i].vec_q, [-1])
    println(vecModelSub[i].mat_d)
    println(vecModelSub[i].vec_q)
    return vecModelSub
end


function setModelSub(mat_a, vec_b, vec_c, vecSense, indexMas, blocks, indexSub)
    numSub = length(blocks)
    vecModelSub = Vector{ModelSub}(undef, numSub)
    for k = 1:numSub
        vecModelSub[k] = ModelSub(
            Model(solver = GurobiSolver(OutputFlag = 0, gurobi_env)),     # mod     !
            deepcopy(mat_a[collect(indexMas), indexSub[k]]),              # mat_e   <- A0
            deepcopy(vec_c[indexSub[k]]),                                 # vec_l
            deepcopy(mat_a[blocks[k], indexSub[k]]),                      # mat_d   <- ASub
            deepcopy(vec_b[blocks[k]]),                                   # vec_q
            0,                                                            # vec_x
            deepcopy(vecSense[blocks[k]])                                 # vecSense
            )
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
