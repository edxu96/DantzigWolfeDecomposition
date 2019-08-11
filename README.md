
# DantzigWolfeDecomposition
Dantzig-Wolfe series of decomposition and reformulation algorithm to solve MILP

This repo has been archived On Aug 11, 2019. New update will be made to [edxu96/MatrixOptim](https://github.com/edxu96/MatrixOptim), which is the aggregation of robust optimization and decomposition.

## 1,  Introduction

According to the PhD thesis by James Richard Tebboth [A Computational Study of Dantzig-Wolfe Decomposition][1]:

> Dantzig-Wolfe decomposition is an optimisation technique for solving large scale, block structured, linear programming (LP) problems. Problems from many different fields such as production planning, refinery optimisation, and resource allocation may be formulated as LP problems. Where there is some structure arising from repeated components in the problem, such as the handling of multi-periods, multi-locations, or multi-products, the problem may potentially be solved using Dantzig-Wolfe decomposition.

> Dantzig-Wolfe decomposition will not rival mainstream techniques as an optimisation method for all LP problems. But we do show that Dantzig-Wolfe decomposition has some niche areas of application: certain large scale classes of primal block angular structured problems, and in particular where the context demands rapid results using parallel optimisation, or near optimal solutions with a guaranteed quality.

> LP optimisation software use two main classes of methods. The simplex method is a gradient descent method that moves along the edge of the feasible region [Chv83, Dan63]. Interior point methods (IPM) move through the interior of the feasible region [Wri97].

DW-decomp is composed of three major parts:
1. Dantzig-Wolfe Reformulation
2. Dantzig-Wolfe Column Generation
3. Stable Column Generation

More detailed explanation can be found in [edxu96/DantzigWolfeDecomposition/wiki](https://github.com/edxu96/DantzigWolfeDecomposition/wiki/1-Home).

## 2,  What can It do Right Now ?   

- [x] Relaxed Mixed Integer Linear Programming Problems (Reformulation)
    - **Input** A linear programming problem with complicating constraints.
    - **Output** The solution of the linear programming problem obtained after using the Dantzig-Wolfe decomposition algorithm.
- [ ] Bounded Mixed Integer Linear Programming (with Extreme Points)
- [ ] Unbounded Mixed Integer Linear Programming (with Extreme Rays)

It must be combined with an efficient branch-and-cut algorithm to solve a MILP.

## 3,  How to use

The code is in five files, which should be stored in your working directory:

```
DantzigWolfeDecomposition.jl
    FuncDW.jl
    FuncSub.jl
    FuncMas.jl
    FuncStab.jl
```

Load other necessary packages:

```Julia
using JuMP
using CPLEX
using Gurobi
using GLPKMathProgInterface
using LinearAlgebra
using MathProgBase
using SparseArrays
```

There is only `doDWDecomp` in the module:

```Julia
using DantzigWolfeDecomposition
## Set parameter and start DW-Decomposition
dualPen = 10
dualPenMult = 0.1
dualPenThreshold = 0.01 - 1e-5
epsilon = 0.00001                # [maximum difference between two bounds]
whePrint = false                 # [whether to print detailed info]
whiSolver = 2                    # 1: Gurobi, 2: CPLEX, 3: GLPK
doDWDecomp(
    mat_a, vec_b, vec_c,                                  # Data in LP Problem
    vecSenseAll, indexMas, blocks, indexSub, numXPerSub,  # Data for DW-Decomp
    dualPen, dualPenMult, dualPenThreshold,               # Para for Stable
    epsilon, whePrint, whiSolver                          # Control Para
    )
```

Attention: The support for `GLPK` is not very well. The `Gurobi` solver does better than `CPLEX`. If you don't have `Gurobi` or `CPLEX`, you may need to comment out the two lines in the `DantzigWolfeDecomposition.jl` file.

```Julia
using CPLEX
using Gurobi
```

In `Test/test.jl`, you can see the example to use `doDWDecomp` to solve Generalized Assignment Problem.

## 4,  Contribution

Edward J. Xu is maintaining the project.

It's originally inspired by Professor Stefan Røpke, DTU Management.

## 5,  Reference

1. Tebboth, J. R., 2001. A computational study of Dantzig-Wolfe decomposition. University of Buckingham.

[1]: http://eaton.math.rpi.edu/CourseMaterials/PreviousSemesters/PreviousSemesters/Spring08/JM6640/tebboth.pdf
