# Dantzig-Wolfe Reformulation and Column Generation
# Edward J. Xu
# 2019.5.28
########################################################################################################################


function getData(whi::Int64)
    if whi == 1  # Convexify􏰁 constraints into one sub-problem
        vec_c = - [292  453  359  219  268  736  291  443  403  498  400  967]
        mat_a = [
            Matrix{Float64}(I, 5, 5) zeros(5, 1) Matrix{Float64}(I, 5, 5) zeros(5, 1);
            14 20 6 16 10 -43 zeros(1, 6);
            zeros(1, 6) 14 20 6 16 10 -43
            ]
        vec_b = hcat([ones(5,1); 1e-6; 1e-6])
        indexSub = [collect(1:12)]
        numSub = 1
        numXInSub = 12
        m_mat_h = 5
        vecRowMatD = [2]
    elseif whi == 2  # Convexify􏰁 constraints into two sub-problems
        vec_c = - [292  453  359  219  268  736  291  443  403  498  400  967]
        mat_a = [
            Matrix{Float64}(I, 5, 5) zeros(5, 1) Matrix{Float64}(I, 5, 5) zeros(5, 1);
            14 20 6 16 10 -43 zeros(1, 6);
            zeros(1, 6) 14 20 6 16 10 -43
            ]
        vec_b = hcat([ones(5,1); 1e-6; 1e-6])
        indexSub = [collect(1:6), collect(7:12)]
        numSub = 2
        numXInSub = 6
        m_mat_h = 5
        vecRowMatD = [1, 1]
    elseif whi == 3  # Convexify􏰁 constraints into one sub-problems, restrain the y
        vec_c = - [292  453  359  219  268  736  291  443  403  498  400  967]
        mat_a = [
            Matrix{Float64}(I, 6, 6) Matrix{Float64}(I, 6, 6);
            14 20 6 16 10 -43 zeros(1, 6);
            zeros(1, 6) 14 20 6 16 10 -43
            ]
        vec_b = hcat([ones(5,1); 2; 1e-6; 1e-6])
        indexSub = [collect(1:12)]
        numSub = 1
        numXInSub = 12
        m_mat_h = 6
        vecRowMatD = [2]
    elseif whi == 4  # Convexify􏰁 constraints into two sub-problems, restrain the y
        vec_c = - [292  453  359  219  268  736  291  443  403  498  400  967]
        mat_a = [
            Matrix{Float64}(I, 6, 6) Matrix{Float64}(I, 6, 6);
            14 20 6 16 10 -43 zeros(1, 6);
            zeros(1, 6) 14 20 6 16 10 -43
            ]
        vec_b = hcat([ones(5,1); 2; 1e-6; 1e-6])
        indexSub = [collect(1:6), collect(7:12)]
        numSub = 2
        numXInSub = 6
        m_mat_h = 6
        vecRowMatD = [1, 1]
    end
    return (vec_c, mat_a, vec_b, indexSub, numSub, numXInSub, m_mat_h, vecRowMatD)
end
