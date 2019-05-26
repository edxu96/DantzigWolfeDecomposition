




function getMat(mod)
    JuMP.build(mod)
    mpb = MathProgBase
    modInternal = internalmodel(mod)
    vecBoundLowerCons = mpb.getconstrLB(modInternal)
    vecBoundUpperCons = mpb.getconstrUB(modInternal)
    vecSense = Vector{Sense}(undef, length(vecBoundLowerCons))
    vec_b = Vector{Float64}(undef, length(vecBoundLowerCons))
    for i = 1:length(vecBoundLowerCons)
        if vecBoundLowerCons[i] == vecBoundUpperCons[i]
            vecSense[i] = eq
            vec_b[i] = vecBoundLowerCons[i]
        elseif vecBoundLowerCons[i] == -Inf
            vecSense[i] = leq
            vec_b[i] = vecBoundUpperCons[i]
        elseif vecBoundUpperCons[i] == Inf
            vecSense[i] = geq
            vec_b[i] = vecBoundLowerCons[i]
        else
            throw("Error: One of the constraints in the input model is a \"range\" constraint. This is not supported.")
        end
    end
    mat_a = mpb.getconstrmatrix(modInternal)
    vec_c = mpb.getobj(modInternal)
    vecBoundLowerVar = mpb.getvarLB(modInternal)
    vecBoundUpperVar = mpb.getvarUB(modInternal)
    vecTypeVar = mpb.getvartype(modInternal)
    return (mat_a, vec_b, vec_c, vecSense)
    # return (mat_a, vec_b, vec_c, vecBoundLowerVar, vecBoundUpperVar, vecTypeVar, vecSense)
end


function addVarsFromRow(Arow, varBitSet)
    for idx=1:length(Arow.nzval)
        if Arow.nzval[idx] != 0
            push!(varBitSet, Arow.nzind[idx])
        end
    end
end


function getIndex(mat_a, blocks)
    # find out which variables that belong to which block:
    numSub = length(blocks)
    if numSub < 1
        throw("error. We expect at least one block")
    end
    (numCons, numVar) = size(mat_a)
    varsInSub = Vector{BitSet}(undef, numSub)
    for k = 1:numSub
        varsInSub[k] = BitSet()
        for row in blocks[k]
            addVarsFromRow(mat_a[row,:], varsInSub[k])
        end
    end
    # check that the blocks are not sharing any variables
    for k1 = 1:numSub
        for k2 = k1 + 1:numSub
            if !isempty(intersect(varsInSub[k1],varsInSub[k2]))
                throw("Blocks $k1 and $k2 share variables. This is not supported.")
            end
        end
    end
    # Check that all variables belong to a block
    setUnion = deepcopy(varsInSub[1])
    for k = 2:numSub
        union!(setUnion, varsInSub[k])
    end
    for j = 1:numVar
        if !(j in setUnion)
            throw("variable $j does not belong to any block. This is not supported")
        end
    end
    # Find the rows that are left in master problem:
    indexMas = BitSet(1:numCons)
    for block in blocks
        setdiff!(indexMas, block)
    end
    ## IndexSub
    indexSub = Vector{Vector{Int64}}(undef, numSub)
    for k = 1:numSub
        indexSub[k] = collect(varsInSub[k])
    end
    return (indexMas, indexSub)
end
