# Dantzig-Wolfe Reformulation and Column Generation
# Functions for Master Problem
# Edward J. Xu
# 2019.5.25
########################################################################################################################


function setModelMas(mat_e, vec_b, num_sub, num_border)
    ## X1: initial matrix of extreme points
    vec_p = vec_b[1: num_border, 1]
    mod_mas = Model(solver = GurobiSolver(OutputFlag = 0, gurobi_env))
    (m_e, n_e) = size(mat_e)
    K = 1
    # In this case we do not use a starting set of extreme points.
    # we just use a dummy starting column
    @variable(mod_mas, vec_lambda[1: K] >= 0 )
    # Remember to consider if we need to maximize or minimize
    @objective(mod_mas, Max, sum(- 1000 * vec_lambda[j] for j = 1: K) )
    # remember to change "==" if your master-problem uses a different type of constraints!
    @constraint(mod_mas, vec_consRef[i = 1: m_e], sum(vec_p[i] * vec_lambda[j] for j = 1: K) == vec_p[i])
    @constraint(mod_mas, vec_consConvex[k = 1: num_sub], sum(vec_lambda[j] for j = 1: K) == 1)
    return (mod_mas, vec_consRef, vec_consConvex, vec_lambda)
end


function solveMas(mod_mas, vec_consRef, vec_consConvex)
    status = solve(mod_mas)
    if status != :Optimal
        throw("Error: Non-optimal master-problem status")
    end
    vec_pi = getdual(vec_consRef)
    vec_kappa = getdual(vec_consConvex)
    obj_master = getobjectivevalue(mod_mas)
    ## ensure that vec_pi and vec_kappa are row vectors
    vec_pi = reshape(vec_pi, 1, length(vec_pi))
    vec_kappa = reshape(vec_kappa, 1, length(vec_kappa))
    return (vec_pi, vec_kappa, obj_master)
end
