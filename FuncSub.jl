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
end


function setVariableSub(modSub::JuMP.Model, mat_d, vec_q)
    (m_d, n_d) = size(mat_d)
    @variable(modSub, vec_x[1: n_d] >= 0, Bin)
    # objective is not here. We define once dual variables become known
    @objective(modSub, Max, 0)
    # remember to change "<=" if your sub-problem uses a different type of constraints!
    @constraint(modSub, cons[i = 1: m_d], sum(mat_d[i, j] * vec_x[j] for j = 1: n_d) <= vec_q[i])
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


function setModelSub(num_sub, num_border, vec_c, mat_a, vec_b)
    vecModelSub = Vector{ModelSub}(undef, num_sub)
    for k = 1: num_sub
        vecModelSub[k] = ModelSub(
            Model(solver = GurobiSolver(OutputFlag = 0, gurobi_env)),     # mod
            mat_a[1: num_border, index_sub[k]],                           # mat_e
            vec_c[index_sub[k]],                                          # vec_l
            hcat(mat_a[(num_border + k), index_sub[k]])',                 # mat_d
            vec_b[(num_border + k), 1],                                   # vec_q
            0                                                             # vec_x
            )
    end
    # Branch
    vecModelSub = setBranch(vecModelSub, 1, 2)
    #
    for k = 1: num_sub
        vecModelSub[k].vec_x = setVariableSub(vecModelSub[k].mod, vecModelSub[k].mat_d, vecModelSub[k].vec_q)
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
    return (costReduce, vec_x_result), 
end
