using Optim
using Interpolations
using Roots
using CSV
using Plots

@time begin
    

println("--------Huggett(1993) 일반균형 추정(bisection) ----------")

# -RHS 값을 시산해주는 함수

function minus_RHS(x, a, e, q, a_grid, V_old, Pi, e_idx, beta)
    sigma = 1.5
    c = a + e - x * q
    v_1 = LinearInterpolation(a_grid, V_old[:,1])
    v_2 = LinearInterpolation(a_grid, V_old[:,2])
    if c <= 0
        RHS_value = 1e4
    else
        RHS_value = -c^(1-sigma) /(1-sigma) - beta * (v_1(x) * Pi[e_idx, 1] + v_2(x)* Pi[e_idx, 2])
    end
    return RHS_value
end;


# a->a' 행렬을 만들어주는 함수

function BTM(a_grid, a_prime)
    # transition matrix 초기화
    n = length(a_grid)
    transition_matrix = zeros(n, n)

    # 각 a_prime 값에 대해 transition matrix 채우기
    for i = 1:length(a_prime)
        # a_prime 값이 a_grid 범위 내에 있는지 확인
        if a_prime[i] >= a_grid[1] && a_prime[i] <= a_grid[end]
            # a_prime이 a_grid 내에 있으면 해당하는 인덱스 찾기
            idx = argmin(abs.(a_grid .- a_prime[i]))

            # transition matrix 업데이트
            if idx <= n # a_prime 값이 a_grid의 범위 내에 있으면
                # a_prime이 a_grid 값과 정확히 일치할 때
                if a_prime[i] == a_grid[idx]
                    transition_matrix[i, idx] = 1
                elseif a_prime[i] > a_grid[idx] # 가까운 쪽이 왼쪽에 있는 경우
                    ratio = (a_prime[i] - a_grid[idx]) / (a_grid[idx+1] - a_grid[idx])
                    transition_matrix[i, idx] = 1 - ratio
                    transition_matrix[i, idx+1] = ratio
                elseif a_prime[i] < a_grid[idx]  # 가까운 쪽이 오른쪽에 있는 경우
                    ratio = (a_prime[i] - a_grid[idx-1]) / (a_grid[idx] - a_grid[idx-1])
                    transition_matrix[i, idx] = ratio
                    transition_matrix[i, idx-1] = 1 - ratio
                end
            else # a_prime 값이 a_grid의 최댓값보다 큰 경우
                transition_matrix[i, i] = 1
            end
        end
    end

    return transition_matrix
end;

function PE_Huggett(q)
    # parameter
    e_grid = [1.0; 0.1]
    beta = 0.99322
    sigma = 1.5
    Pi = [0.925    1-0.925;
          1-0.5    0.5]
    a_lb = -2
    a_ub = 4
    N = 150    # number of a_grid
    K = 2      # number of e_grid
    a_grid = range(a_lb, stop=a_ub, length=N);
    
    println("-------------------------------")
    println("채권가격 q: ", q)
    
    ## (Step1) VFI로부터 policy function 얻기
    
    # prepare for the PFI
    V_ini = zeros(N, K);             # intial guess for the policy function
    V_old = zeros(N, K);
    V_new = zeros(N, K);
    pol_a_prime = zeros(N, K);
    V_old .= V_ini
    
    # loop
    max_iter = 5000
    count = 0
    tol = 1e-4
    diff = 10
    
    while count < max_iter && diff > tol
        for e_idx=1:K
            for a_idx=1:N
                a = a_grid[a_idx]
                e = e_grid[e_idx]
                function f_obj(x)
                    obj_value = minus_RHS(x, a, e, q, a_grid, V_old, Pi, e_idx, beta)
                    return obj_value
                end
                # a_lb <= a' <= (a+e)/q (C가 음수가 안될 조건)
                a_ub_endo = min((a+e)/q, a_ub)
                results = optimize(f_obj, a_lb, a_ub_endo)
                pol_a_prime[a_idx, e_idx] = Optim.minimizer(results)
                V_new[a_idx, e_idx] = -Optim.minimum(results)
                if Optim.converged(results) && Optim.minimum(results) == 1e4
                    pol_a_prime[a_idx, e_idx] = a_lb
                end
            end
        end
        diff = maximum(abs.(V_new .- V_old))
        count += 1
        copyto!(V_old, V_new)
    end
    
    # (Step2) Stationary distribution 얻기
    
    # psi_0 (initial guess for the stationary distribution)
    psi_ini = zeros(N*K, 1)
    psi_ini[1] = 1
    
    a_h_tr = BTM(a_grid, pol_a_prime[:, 1])
    a_l_tr = BTM(a_grid, pol_a_prime[:, 2])
    W = [Pi[1, 1]*a_h_tr Pi[1, 2]*a_h_tr;
         Pi[2, 1]*a_l_tr Pi[2, 2]*a_l_tr]
    
    step2_count = 0
    step2_diff = 10
    step2_tol = 1e-6
    psi_old = copy(psi_ini)
    psi_new = zeros(N*K, 1)
    
    while step2_count < max_iter && step2_diff > step2_tol
        psi_new .= W' * psi_old  # Notation 주의!  
        step2_diff = maximum(abs.(psi_new .- psi_old))
        step2_count += 1
        psi_old .= psi_new
    end
    
    # CDF 계산
    psi_mat = hcat(psi_new[1:N], psi_new[N+1:end])
    Psi_mat = copy(psi_mat)
    for i=2:N
        Psi_mat[i, :] .= Psi_mat[i-1, :] .+ psi_mat[i, :]
    end
    
    ## (Step 3) Update
    
    # aggregation
    A = sum(reshape(a_grid,1,:)*psi_mat)
    
    println("총 채권 수요(A) : ", A) 
 
    return A, pol_a_prime, psi_mat, Psi_mat
end

q_sol = fzero(PE_Huggett, (0.98, 1.02), Roots.Brent(), xtol=1e-6)
println("균형 채권가격(q*): ", q_sol )

A, pol_a_prime, psi_mat, Psi_mat = PE_Huggett(q_sol)

# 최종 결과 저장
# 행렬을 2차원 튜플의 배열로 변환
pol_a_sol = [Tuple(row) for row in eachrow(pol_a_prime)]
density_sol = [Tuple(row) for row in eachrow(psi_mat)]
distribution_sol = [Tuple(row) for row in eachrow(Psi_mat)]

# 파일 경로
file_path1 = "C:\\Users\\82106\\Documents\\4. 2024학년도\\Road to HANK\\(2) Julia\\policy.csv"
file_path2 = "C:\\Users\\82106\\Documents\\4. 2024학년도\\Road to HANK\\(2) Julia\\density.csv"
file_path3 = "C:\\Users\\82106\\Documents\\4. 2024학년도\\Road to HANK\\(2) Julia\\distribution.csv"

# 데이터를 CSV 파일에 쓰기
CSV.write(file_path1, pol_a_sol, header=false)
CSV.write(file_path2, density_sol, header=false)
CSV.write(file_path3, distribution_sol, header=false)


# 그림 생성

# parameter
a_lb = -2
a_ub = 4
N = 150    # number of a_grid
K = 2      # number of e_grid
   a_grid = range(a_lb, stop=a_ub, length=N);

# Plot 생성
# 그림1: ploicy function
p1 = plot(a_grid, pol_a_prime[:, 1], label="e=1", linecolor=:blue, xlabel="a_grid", ylabel="a(a,e)", lw=2)
plot!(p1, a_grid, pol_a_prime[:, 2], label="e=0.1", linecolor=:red, lw=2)
plot!(p1, a_grid, a_grid, label="diagonal", linecolor=:green, linestyle=:dash, lw=2)
savefig(p1, "C:\\Users\\82106\\Documents\\4. 2024학년도\\Road to HANK\\(2) Julia\\policy.png") 

# 그림2: distribution
p2 = plot(a_grid, Psi_mat[:, 1], label="e=1", linecolor=:blue, xlabel="a_grid", ylabel="distribution", lw=2)
plot!(p2, a_grid, Psi_mat[:, 2], label="e=0.1", linecolor=:red, lw=2)
savefig(p2, "C:\\Users\\82106\\Documents\\4. 2024학년도\\Road to HANK\\(2) Julia\\distribution.png") 

# 그림3: density
p3 = plot(a_grid, psi_mat[:, 1], label="e=1", linecolor=:blue, xlabel="a_grid", ylabel="density", lw=2)
plot!(p3, a_grid, psi_mat[:, 2], label="e=0.1", linecolor=:red, lw=2)
savefig(p3, "C:\\Users\\82106\\Documents\\4. 2024학년도\\Road to HANK\\(2) Julia\\density.png") 

end

