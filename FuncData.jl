

function getStrNameVar(numSub, numXPerSub)
    vecStrNameVar = Vector{String}(undef, numSub * numXPerSub)
    idx = 1
    for i = 1:numSub
        for j = 1:numXPerSub
            vecStrNameVar[idx] = "x_($i, $j)"
            idx += 1
        end
    end
    return vecStrNameVar
end


function setGAPModel(w, p, c)
    (numSub, numXPerSub) = size(w)
    # mod = Model(solver = CplexSolver(CPX_PARAM_TILIM = 1))
    modAll = Model(solver = GurobiSolver())
    @variable(modAll, x[1:numSub, 1:numXPerSub], Bin)
    @objective(modAll, Min, sum(p[i,j] * x[i,j] for i = 1:numSub for j = 1:numXPerSub))
    # each job must be served
    @constraint(modAll, [j = 1:numXPerSub], sum(x[i,j] for i = 1:numSub) == 1)
    @constraint(modAll, [i = 1:numSub], sum(w[i,j] * x[i,j] for j = 1:numXPerSub) <= c[i])
    ##
    # define blocks (Each block becomes a sub-problem)
    blocks = []
    for i = 1:numSub
        # constraint (numXPerSub + i) goes into block i
        push!(blocks, [numXPerSub + i])
    end
    # When using cplex we need to start solving the model before we can
    # export the constraint matrix and so on.
    # solve(mod)
    JuMP.build(modAll)
    return (modAll, numXPerSub, blocks)
end


function readData(strFileName::String, series::Int64)
    # function for reading GAP files from the OR-library
    open(strFileName) do f
        allFile = read(f, String)
        allNumbers = map(x -> parse(Int64,x), split(allFile))
        # return allNumbers
        numSeries = allNumbers[1]
        println("File contains $numSeries instances.")
        if series > numSeries
            throw("readData(...): specified file does not contain $(series) example")
        end
        index = 2
        for inst = 1:numSeries
            numSub = allNumbers[index]
            numXPerSub = allNumbers[index + 1]
            println("numSub = $(numSub), numXPerSub = $(numXPerSub)")
            w = zeros(numSub, numXPerSub)
            p = zeros(numSub, numXPerSub)
            c = zeros(numSub)
            index += 2
            println("start index when reading p: $index")
            for i = 1:numSub
                for j = 1:numXPerSub
                    p[i,j] = allNumbers[index]
                    index += 1
                end
            end
            println("start index when reading w: $index")
            for i = 1:numSub
                for j = 1:numXPerSub
                    w[i,j] = allNumbers[index]
                    index += 1
                end
            end
            println("start index when reading c: $index")
            for i=1:numSub
                c[i] = allNumbers[index]
                index += 1
            end
            if inst == series
                return (w, p, c)
            end
        end
    end
end
