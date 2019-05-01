# Dantzig-Wolfe Reformulation and Column Generation
# Edward J. Xu
# 2019.5.1
push!(LOAD_PATH, "$(homedir())/Desktop")
cd("$(homedir())/Desktop")
using JuMP
# using CPLEX
using Gurobi


mutable struct ModelSub
    mod::JuMP.Model
    mat_e
    vec_l
    mat_d
    vec_q
    vec_x
end


function setVariableSub(mod_sub::JuMP.Model, mat_d, vec_q)
    (m_d, n_d) = size(mat_d)
    @variable(mod_sub, vec_x[1: n_d] >= 0, Bin)
    # objective is not here. We define once dual variables become known
    @objective(mod_sub, Max, 0)
    # remember to change "<=" if your sub-problem uses a different type of constraints!
    @constraint(mod_sub, cons[i = 1: m_d], sum(mat_d[i, j] * vec_x[j] for j = 1: n_d) <= vec_q[i])
    return vec_x
end


function setModelSub(num_sub, vec_c, mat_a, vec_b)
    modelSub = Vector{ModelSub}(undef, num_sub)
    for k = 1: num_sub
        modelSub[k] = ModelSub(
            Model(solver = GurobiSolver(OutputFlag = 0, gurobi_env)),   # mod
            mat_a[1: 3, index_sub[k]],                                  # mat_e
            vec_c[index_sub[k]],                                        # vec_l
            hcat(mat_a[(3 + k), index_sub[k]])',                        # mat_d
            vec_b[(3 + k), 1],                                          # vec_q
            0
            )
    end
    for k = 1: num_sub
        modelSub[k].vec_x = setVariableSub(modelSub[k].mod, modelSub[k].mat_d, modelSub[k].vec_q)
    end
    return modelSub
end


# X1: initial matrix of extreme points
function setModelMas(mat_e, vec_b, num_sub)
    vec_p = vec_b[1: 3, 1]
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
    return (mod_mas, vec_consRef, vec_consConvex)
end


function addColToMaster(mod_mas::JuMP.Model, mod_sub, vec_x_result, vec_consRef, consConvex)
    (m, n) = size(mod_sub.mat_e)
    A0x = mod_sub.mat_e * vec_x_result
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
    objCoef = sum(mod_sub.vec_l[j] * vec_x_result[j] for j = 1: n)
    println("objCoef = $objCoef.")
    @variable(
        mod_mas,
        lambdaNew >= 0,                     # New variable to be added
        objective = objCoef,                # cost coefficient of new varaible in the objective
        inconstraints = vec_consTouched,    # constraints to be modified
        coefficients = vals                 # the coefficients of the variable in those constraints
        )
end


function solveSub(mod_sub, vec_pi, kappa)
    (m, n) = size(mod_sub.mat_e)
    piA0 = vec_pi * mod_sub.mat_e
    @objective(mod_sub.mod, Max, sum(mod_sub.vec_l[j] * mod_sub.vec_x[j] for j = 1: n) -
        sum(piA0[1, j] * mod_sub.vec_x[j] for j = 1: n) - kappa)
    status = solve(mod_sub.mod)
    if status != :Optimal
        throw("Error: Non-optimal sub-problem status")
    end
    return (getobjectivevalue(mod_sub.mod), getvalue(mod_sub.vec_x))
end


function main()
    timeStart = time()
    epsilon = 0.00001
    global gurobi_env = Gurobi.Env()
    ## 1,  Data
    include("Data_GAP_1.jl")
    ## 2,  Set Up Sub Models
    modelSub = setModelSub(num_sub, vec_c, mat_a, vec_b)
    ## 3,  Set Up Master Models
    (mod_mas, vec_consRef, vec_consConvex) = setModelMas(modelSub[1].mat_e, vec_b, num_sub)
    ## 4,  Begin Iteration
    done = false
    iter = 1
    while !done
        status = solve(mod_mas)
        if status != :Optimal
            throw("Error: Non-optimal master-problem status")
        end
        vec_pi = getdual(vec_consRef)
        # ensure that vec_pi and vec_kappa are  row vectors
        vec_pi = reshape(vec_pi, 1, length(vec_pi))
        vec_kappa = getdual(vec_consConvex)
        vec_kappa = reshape(vec_kappa, 1, length(vec_kappa))
        println("vec_pi = $(vec_pi)\n",
                "vec_kappa = $(vec_kappa)")
        done = true
        #
        println("--------------------------------------------------------------------------------\n",
                "$(iter)-th iteration. obj_master = $(getobjectivevalue(mod_mas)).\n",
                "--------------------------------------------------------------------------------")
        costReduceBest = -1
        for k = 1: num_sub
            (costReduce, vec_x_result) = solveSub(modelSub[k], vec_pi, vec_kappa[k])
            println("Reduced cost of $(k)-th sub model is $costReduce.")
            if costReduce > costReduceBest
                costReduceBest = costReduce
            end
            # remember to look for negative reduced cost if we are minimizing.
            if costReduce > epsilon
                addColToMaster(mod_mas, modelSub[k], vec_x_result, vec_consRef, vec_consConvex[k])
                done = false
            end
        end
        iter += 1
        println("best reduced cost is $costReduceBest")
    end
    ## 5,  Print the Result
    println("################################################################################\n",
            "Done after $(iter-1) iterations.\n",
            "obj_master = $(getobjectivevalue(mod_mas))\n",
            "################################################################################")
    println("Elapsed time is $(time() - timeStart) seconds.")
end


main()
