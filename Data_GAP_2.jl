# Dantzig-Wolfe Reformulation and Column Generation
# Data for General Assignment Problem 1
# Edward J. Xu
# 2019.5.1


function getData()
    vec_c = [#=
        =# 9 9 1 4 6 4 2 3 5 1 9 9 7 6 3  #=
        =# 9 1 8 6 7 8 2 6 6 5 3 4 7 5 3  #=
        =# 4 8 1 8 1 4 1 4 2 6 2 1 2 5 1  #=
        =# 7 9 7 8 3 6 2 3 4 7 9 9 3 9 1]
    mat_a = [
        Matrix{Float64}(I, 15, 15) Matrix{Float64}(I, 15, 15) Matrix{Float64}(I, 15, 15) Matrix{Float64}(I, 15, 15);
        8 6 1 7 7 7 5 5 7 7 3 3 8 5 4 zeros(1,45) ;
        zeros(1,15) 5 5 3 8 5 6 5 9 9 6 1 9 5 6 3 zeros(1,30);
        zeros(1,30) 3 2 4 4 9 1 7 3 3 3 5 3 7 7 1 zeros(1,15);
        zeros(1,45) 6 7 9 9 3 5 2 1 5 5 4 4 6 2 1
        ]
    vec_b = hcat([ones(15,1); 22.0; 18.0; 18.0; 19.0])
    index_sub = [collect(1: 15), collect(16: 30), collect(31: 45), collect(46: 60)]
    num_sub = 4
    num_x_sub = 15
    num_border = 15
    #
    vecStr_nameVar = Vector{String}(undef, length(vec_c))
    idx = 1
    for i = 1: num_sub
        for j = 1: num_x_sub
            vecStr_nameVar[idx] = "x_{$(i),$(j)}"
            idx += 1
        end
    end
    return (vec_c, mat_a, vec_b, index_sub, num_sub, num_x_sub, num_border, vecStr_nameVar)
end


(vec_c, mat_a, vec_b, index_sub, num_sub, num_x_sub, num_border, vecStr_nameVar) = getData()
