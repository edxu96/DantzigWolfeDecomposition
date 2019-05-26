# DantzigWolfeDecomposition

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

## 2,  How to use

The code is in four files:

```
FuncDW.jl
    FuncSub.jl
    FuncMas.jl
    FuncStab.jl
```

To use the function `doDWDecomp`:

```Julia
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
```

In `main.jl`, you can see the example to use `doDWDecomp` in solving Generalized Assignment Problem.

## 3,  Contributions

Edward J. Xu is maintaining the project.

It's originally inspired by Professor Stefan Røpke, DTU Management.
