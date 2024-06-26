# Housing model with a Life cycle 
# Yeongwoong Do (2024-04-26)

using LinearAlgebra, Interpolations, Plots, Printf, Dates, MAT, Optim, SparseArrays, NLopt, Roots, DelimitedFiles, Statistics, NLsolve, ProgressMeter

###############################################################
###################### (subfunctions) #########################
###############################################################

# income tax 
function income_tax(y)
    # income tax rate
    # mean earning = 43,330,000 
    t_y     = [0.06, 0.15, 0.24, 0.35, 0.38]
    thres_y = [0.52, 1.45, 2.37, 3.83, 7.37]
    d_y     = [0.3, 0.6, 0.85, 0.95, 0.98]
    d_thres = [0.12, 0.35, 1.04, 2.31]
    tax     = [0, 0.0194, 0.1440, 0.3545, 0.8553]
    if y <= thres_y[1]
        if y<= d_thres[1]
            return (y*d_y[1])*t_y[1]
        elseif y<= d_thres[2]
            return (d_thres[1]*d_y[1] + (y-d_thres[1])*d_y[2])*t_y[1]
        else
            return (d_thres[1]*d_y[1] + (d_thres[2]-d_thres[1])*d_y[2]+(y-d_thres[2])*d_y[3])*t_y[1]
        end
    elseif y <= thres_y[2]
        if y<= d_thres[3]
            return tax[2] + ((y-thres_y[1])*d_y[3])*t_y[2]
        else
            return tax[2] + ((d_thres[3]-thres_y[1])*d_y[3]+(y-d_thres[3])*d_y[4])*t_y[2]
        end
    elseif y <= thres_y[3]
        if y<= d_thres[4]
            return tax[3] + ((y-thres_y[2])*d_y[4])*t_y[3]
        else
            return tax[3] + ((d_thres[4]-thres_y[2])*d_y[4]+(y-d_thres[4])*d_y[5])*t_y[3]
        end
    elseif y <= thres_y[4]
        return tax[4] + (y-thres_y[3])*d_y[5]*t_y[4]
    else
        return tax[5] + (y-thres_y[4])*d_y[5]*t_y[5]
    end
end;

# housing tax
function housing_tax_cost(gamma, delta, h, h_prime, cutoff, q, tau_p, HRET)
    transaction_cost = gamma*(q*h_prime + q*(1-delta)*h)*(h_prime!=h)
    property_tax = tau_p*q*h_prime
    HRET_tax = HRET*(q*h_prime-cutoff)*(q*h_prime>=cutoff)
    return transaction_cost + property_tax + HRET_tax
end

# RHS of Renter's problem
function RHS_renter(a_prime, bar_x, sigma, alpha, beta, V_itp, Pi, s_prob, rho, z_idx, t, b_bar, gamma, q, h, delta) 
    x = bar_x - a_prime -gamma*q*(1-delta)*h
    if x <0
        return -1e10
    end
    u = (x^(1-sigma)/(1-sigma))*(1/rho)^((1-alpha)*(1-sigma))
    T = length(s_prob)
    if t==T   
        return u + b_bar*(a_prime)^(1-sigma)/(1-sigma)
    else
        exp = 0
        K0, K = size(Pi)
        @inbounds for k=1:K
            # 딕셔너리에서 보간 함수를 가져옴
            itp = V_itp["func1_$(k)"]
            exp = exp + Pi[z_idx,k]*itp(a_prime)
        end 
        return u + beta*s_prob[t]*exp
    end
end;

# RHS of onwer-occupier's problem
function RHS_owner(a_prime, h_prime, h_prime_idx, bar_x, sigma, alpha, beta, V_itp, Pi, s_prob, q, z_idx, t, h_cost, b_bar) 
    c = bar_x -a_prime - h_cost - q*h_prime
    if c <=0 || h_prime<=0
        return -1e10
    end
    T = length(s_prob)
    u = ((c/alpha)^alpha * (h_prime/(1-alpha))^(1-alpha))^(1-sigma)/(1-sigma)
    if t==T
        return u + b_bar*(a_prime+q*h_prime)^(1-sigma)/(1-sigma)
    else
        exp = 0; K0, K = size(Pi);
        @inbounds for k=1:K
            itp = V_itp["func$(h_prime_idx)_$(k)"]
            exp = exp + Pi[z_idx,k]*itp(a_prime)
        end       
        return u + beta*s_prob[t]*exp
    end
end;

# RHS of landlord's problem
function RHS_landlord(a_prime, h_prime, h_prime_idx, bar_x, sigma, alpha, beta, V_itp, Pi, s_prob, rho, q, z_idx, t, tau_R, h_cost, b_bar, phi, phi_0)
    x = bar_x -a_prime - h_cost - (q - rho*(1-tau_R)+ phi)*h_prime - phi_0
    s = (1-alpha)*x/(rho*(1-tau_R)-phi)
    c = s*(rho*(1-tau_R)-phi)*alpha/(1-alpha)
    if c<0 || s<0        # c> 0, s>0
        return -1e10
    end
    if h_prime <= s  # h'>s
        return -1e10
    end
    T = length(s_prob)
    u = (x^(1-sigma)/(1-sigma))*(1/(rho*(1-tau_R)))^((1-alpha)*(1-sigma))
    if t==T
        return u + b_bar*(a_prime+q*h_prime)^(1-sigma)/(1-sigma)
    else
        exp = 0; K0, K = size(Pi);
        @inbounds for k=1:K
            itp = V_itp["func$(h_prime_idx)_$(k)"]
            exp = exp + Pi[z_idx,k]*itp(a_prime)
        end       
        return u + beta*s_prob[t]*exp
    end
end;

# aggregation을 위한 transtion matrix 생성 함수
function BTM(pol_a, pol_h, Pi, N,J,K, a_grid, h_grid)
    # 1단계 a, h의 transtion block 만들기
    block_transition = zeros(N*J,N*J,K)
    for k=1:K
        for n=1:N
            for j=1:J
            
                # 값 알아보기 쉽게하기
                h_prime = pol_h[n,j,k]
                a_prime = pol_a[n,j,k]
            
                # h->h' 위치찾기
                h_idx = argmin(abs.(h_grid .- h_prime))                
                # a->a' 위치찾기
                a_idx = argmin(abs.(a_grid .- a_prime))           
                row_num = (j-1)*N + n
                col_num = (h_idx-1)*N + a_idx
                
                if a_prime == a_grid[a_idx] # a'이 a 그리드값과 동일
                    block_transition[row_num, col_num, k] = 1
                
                elseif a_prime > a_grid[a_idx] # a그리드랑 가까운 쪽이 왼쪽인 경우
                    ratio = (a_prime - a_grid[a_idx]) / (a_grid[a_idx+1] - a_grid[a_idx])
                    block_transition[row_num, col_num, k] = 1 - ratio
                    block_transition[row_num, col_num+1, k] = ratio
                
                elseif a_prime < a_grid[a_idx] # a그리드랑 가까운 쪽이 오른쪽인 경우
                    ratio = (a_grid[a_idx] - a_prime) / (a_grid[a_idx] - a_grid[a_idx-1])
                    block_transition[row_num, col_num, k] = 1- ratio
                    block_transition[row_num, col_num-1, k] = ratio
                end
            
            end
        end    
    end
    
    
    # 2단계 block matrix를 배치하기
    W = zeros(N*J*K,N*J*K);

    for k1=1:K
        for k2=1:K
            # block indexing
            r_start = N*J*(k1-1)+1
            r_end   = N*J*k1
            c_start = N*J*(k2-1)+1
            c_end   = N*J*k2
            W[r_start:r_end, c_start:c_end] = Pi[k1,k2].*block_transition[:,:,k1];
        end
    end
    
    return W
end;

# Gini 계산기
function Gini(pop,x, tol)
    # pop의 원소가 0보다 큰 인덱스만 필터링
    # 이것을 하지 않으면 아무도 비싼 자산을 안가지고 있는데, 가지고 있는 것처럼 나와서 지니계수가 이상하게 나옴
    valid_idx = findall(p -> p > tol, pop)
    # 유효한 인덱스를 사용하여 pop과 x 벡터를 필터링
    filtered_pop = pop[valid_idx]
    filtered_x = x[valid_idx]

    # x를 기준으로 정렬된 인덱스를 얻기
    idx = sortperm(filtered_x)
    # 인덱스를 사용하여 두 벡터 모두 정렬
    sorted_x = filtered_x[idx]
    sorted_pop = filtered_pop[idx]
    
    # gini 계수 시산
    N = length(filtered_pop)
    sorted_x = sorted_x./sum(sorted_x)
    sorted_pop = sorted_pop./sum(sorted_pop)
    cum_pop = zeros(N)
    cum_x = zeros(N)
    cum_pop[1] = sorted_pop[1]
    cum_x[1] = sorted_x[1]
    under_area = sorted_x[1]*sorted_pop[1]*0.5
    for i=2:N
        cum_pop[i]=cum_pop[i-1] + sorted_pop[i]
        cum_x[i]=cum_x[i-1] + sorted_x[i]
        under_area = under_area + 0.5*sorted_pop[i]*sorted_x[i] + cum_x[i-1]*sorted_pop[i]
    end
    gini = (0.5-under_area)/0.5
    return gini
end;


#############################################################
####################(Main)###################################
#############################################################

function Housing_PE_LC(q, rho, Tr, phi, phi_0, gamma, delta, h_grid, cutoff, b_bar, HRET, Tr_r, ife_mat, Idio_mat)
    println("-----------------------------------------------------------------------------")
   
    # save vector
    w_age  = Life_mat[:,1]
    T = length(w_age)
    s_prob = Life_mat[:,2]
    s_prob[end] = 0
    mu     = Life_mat[:,3]
    z_grid = Idio_mat[:,1]
    K = length(z_grid)
    Pi     = Idio_mat[:,2:K+1]
    stat_dist = Idio_mat[:,end]

    # parameters

    sigma   = 2.000;       # inverse of EIS
    alpha   = 0.778;       # share of non-housing consumption
    r       = 0.024;       # interest rate
    kappa   = 0.013;       # spread of mortgage debt
    beta    = 0.975;       # time discount factor
    LTV     = 0.600;       # LTV limit
    tau_p   = 0.001;       # property tax rate (effective rate)
    pen     = 0.149;       # pension 
    R       = 37;       # retire age
    n       = 0;
    minus_inf = -1e10;

    # A grid, H_grid
    a_min = -10
    a_max = 10
    N = 100
    J = length(h_grid)
    # a_grid 생성
    a_debt_grid = -(LinRange(sqrt(-a_min), sqrt(0), Int(N/2))).^2
    a_delta = a_debt_grid[end] - a_debt_grid[end-1]
    a_saving_grid = exp.(LinRange(log(a_delta), log(a_max), Int(N/2)))
    a_grid = [a_debt_grid; a_saving_grid]
    a_lb_renter = 0

    # full grid
    A_grid = zeros(N,J,K,T)
    H_grid = zeros(N,J,K,T)
    Y_grid = zeros(N,J,K,T)

    # 빈공간 만들기
    V = zeros(N,J,K,T);
    pol_stance = zeros(N,J,K,T);
    pol_a = zeros(N,J,K,T);
    pol_h = zeros(N,J,K,T);
    pol_s = zeros(N,J,K,T);
    pol_c = zeros(N,J,K,T);
    V_itp = Dict()
    Threads.@threads for k=1:K
        Threads.@threads for j=1:J
            interp = LinearInterpolation(a_grid, V[:,j,k,end])
            V_itp["func$(j)_$(k)"] = interp
        end
    end

@time begin
    
    @showprogress for t=T:-1:1  # Backward!  
            
            Threads.@threads for z_idx=1:K
                Threads.@threads for h_idx=1:J
                    Threads.@threads for a_idx=1:N
                        
                        z = z_grid[z_idx]
                        h = h_grid[h_idx]
                        H_grid[a_idx, h_idx, z_idx, t] = h_grid[h_idx]
                        a = a_grid[a_idx]
                        A_grid[a_idx, h_idx, z_idx, t] = a_grid[a_idx]
                        y = w_age[t]*z
                        Y_grid[a_idx, h_idx, z_idx, t] = w_age[t]*z
                        tax_base = y + r*a*(a>0)

                        # initializing the W()
                        W_R = 0
                        W_R_pol_a = 0
                        W_R_pol_s = 0
                        W_R_pol_c = 0
                        W_o = 0
                        W_o_pol_a = 0
                        W_o_pol_h = 0
                        W_o_pol_s = 0
                        W_o_pol_c = 0
                        W_l = 0
                        W_l_pol_a = 0
                        W_l_pol_h = 0
                        W_l_pol_s = 0
                        W_l_pol_c = 0

                        bar_x = y + (1+r)*a  + kappa*a*(a<0) + q*(1-delta)*h -income_tax(tax_base) + pen*(t>=R) + Tr  # maximum expenditure
                        inc_tax_rate = income_tax(tax_base)/tax_base

                        ###########################################
                        # Renter's problem W(.,h'=0)
                        
                        a_ub_renter = min(bar_x-gamma*q*(1-delta)*h+Tr_r , a_max)
                        if a_ub_renter < a_lb_renter # No solution case
                            W_R = minus_inf
                            W_R_pol_a = 0
                            W_R_pol_s = 0 
                            W_R_pol_c = 0 
                        else
                            function f_rent(x)
                                return -RHS_renter(x, bar_x+Tr_r, sigma, alpha, beta, V_itp, Pi, s_prob, rho, z_idx, t, b_bar, gamma, q, h, delta) 
                            end
                            results = Optim.optimize(f_rent, a_lb_renter, a_ub_renter)
                            W_R = -Optim.minimum(results)
                            W_R_pol_a = Optim.minimizer(results)
                            W_R_pol_s = (1-alpha)*(bar_x+Tr_r- W_R_pol_a)/rho
                            W_R_pol_c = alpha*(bar_x+Tr_r- W_R_pol_a)
                        end
                        ################################################
                        
                        ############################################
                        # Owner-occupier's problem W(.,h'=s)
                        W_o_temp = zeros(J-1);
                        pol_o_a_temp = zeros(J-1);
            
                        Threads.@threads for h_prime_idx = 2:J
                            h_prime = h_grid[h_prime_idx]
                            h_cost = housing_tax_cost(gamma, delta, h, h_prime, cutoff, q, tau_p, HRET)
                            if t == T
                                a_lb_owner = 0  # t==T일 경우, 0을 사용
                            else
                                a_lb_owner = max(-LTV*q*h_prime, a_min)
                            end
                            a_ub_owner = min(bar_x - q*h_prime - h_cost, a_max)

                            if a_ub_owner < a_lb_owner # 말도 안되는 상황이 오면
                                W_o_temp[h_prime_idx-1]= minus_inf
                                pol_o_a_temp[h_prime_idx-1]= 0
                            else
                                function f_owner(x)
                                    return - RHS_owner(x, h_prime, h_prime_idx, bar_x, sigma, alpha, beta, V_itp, Pi, s_prob, q, z_idx, t, h_cost, b_bar)
                                end
                                results_owner_occupied = Optim.optimize(f_owner, a_lb_owner, a_ub_owner)
                                W_o_temp[h_prime_idx-1] = -Optim.minimum(results_owner_occupied)
                                pol_o_a_temp[h_prime_idx-1] = Optim.minimizer(results_owner_occupied)
                            end
                        end

                        W_o = maximum(W_o_temp)
                        o_idx = argmax(W_o_temp)
                        W_o_pol_a = pol_o_a_temp[o_idx]
                        W_o_pol_h = h_grid[o_idx+1]
                        W_o_pol_s = W_o_pol_h
                        W_o_pol_c = bar_x - W_o_pol_a - q*W_o_pol_h - housing_tax_cost(gamma, delta, h, W_o_pol_h, cutoff, q, tau_p, HRET)

                        ################################################################

                        ################################################################
                        # Landlord's problem W(.,h'>s)
                        W_l_temp = zeros(J-1);
                        pol_l_a_temp = zeros(J-1);
                        pol_l_s_temp = zeros(J-1);   

                        Threads.@threads for h_prime_idx = 2:J
                            h_prime = h_grid[h_prime_idx]
                            h_cost = housing_tax_cost(gamma, delta, h, h_prime, cutoff, q, tau_p, HRET)
                            if t == T
                                a_lb_landlord = 0  # t==T일 경우, 0을 사용
                            else
                                a_lb_landlord = max(-LTV*q*h_prime, a_min)
                            end
                            a_ub_landlord = min(bar_x - q*h_prime - h_cost - phi_0, a_max)

                            if a_ub_landlord < a_lb_landlord # 말도 안되는 상황이 오면
                                W_l_temp[h_prime_idx-1]=minus_inf
                                pol_l_a_temp[h_prime_idx-1]= 0
                                pol_l_s_temp[h_prime_idx-1]= 0
                            else
                                function f_landlord(x)
                                    return - RHS_landlord(x, h_prime, h_prime_idx, bar_x, sigma, alpha, beta, V_itp, Pi, s_prob, rho, q, z_idx, t, inc_tax_rate, h_cost, b_bar, phi, phi_0)
                                end
                                results_landlord = Optim.optimize(f_landlord, a_lb_landlord, a_ub_landlord)
                                W_l_temp[h_prime_idx-1] = -Optim.minimum(results_landlord)
                                pol_l_a_temp[h_prime_idx-1] = Optim.minimizer(results_landlord)
                                pol_l_s_temp[h_prime_idx-1] = (1-alpha)*(bar_x -pol_l_a_temp[h_prime_idx-1] - h_cost - (q-rho*(1-inc_tax_rate)+phi)*h_prime- phi_0)/(rho*(1-inc_tax_rate)-phi)
                            end

                        end

                        W_l = maximum(W_l_temp)
                        l_idx = argmax(W_l_temp)
                        W_l_pol_a = pol_l_a_temp[l_idx]
                        W_l_pol_h = h_grid[l_idx+1]
                        W_l_pol_s = pol_l_s_temp[l_idx]
                        W_l_pol_c = W_l_pol_s*(rho*(1-inc_tax_rate)-phi)*alpha/(1-alpha)

                        if W_l_pol_h <= W_l_pol_s && W_l!= minus_inf
                            println("error ! landlord contradiction, H_prime :, $W_l_pol_h, S: $W_l_pol_s, W_l: $W_l")
                        end

                        ###################################        
                        
                        if W_R == max(W_R, W_o, W_l) # HH will be a renter
                            V[a_idx, h_idx, z_idx, t]=W_R
                            pol_stance[a_idx, h_idx, z_idx, t]=0
                            pol_a[a_idx, h_idx, z_idx, t] = W_R_pol_a
                            pol_h[a_idx, h_idx, z_idx, t] = 0
                            pol_s[a_idx, h_idx, z_idx, t] = W_R_pol_s
                            pol_c[a_idx, h_idx, z_idx, t] = W_R_pol_c
                        elseif W_o ==  max(W_R, W_o, W_l) # HH will be a owner-occupier
                            V[a_idx, h_idx, z_idx, t]=W_o
                            pol_stance[a_idx, h_idx, z_idx, t]=1
                            pol_a[a_idx, h_idx, z_idx, t] = W_o_pol_a
                            pol_h[a_idx, h_idx, z_idx, t] = W_o_pol_h
                            pol_s[a_idx, h_idx, z_idx, t] = W_o_pol_s
                            pol_c[a_idx, h_idx, z_idx, t] = W_o_pol_c
                        elseif W_l == max(W_R, W_o, W_l) # HH will be a landlord
                            V[a_idx, h_idx, z_idx, t]=W_l
                            pol_stance[a_idx, h_idx, z_idx, t]=2
                            pol_a[a_idx, h_idx, z_idx, t] = W_l_pol_a
                            pol_h[a_idx, h_idx, z_idx, t] = W_l_pol_h
                            pol_s[a_idx, h_idx, z_idx, t] = W_l_pol_s
                            pol_c[a_idx, h_idx, z_idx, t] = W_l_pol_c
                        elseif W_R==W_o && W_o == W_l   # indifferent = renter
                            V[a_idx, h_idx, z_idx, t]=W_R
                            pol_stance[a_idx, h_idx, z_idx, t]=0
                            pol_a[a_idx, h_idx, z_idx, t] = W_R_pol_a
                            pol_h[a_idx, h_idx, z_idx, t] = 0
                            pol_s[a_idx, h_idx, z_idx, t] = W_R_pol_s
                            pol_c[a_idx, h_idx, z_idx, t] = W_R_pol_c
                        else
                            println("error !  W_R: $W_R, W_o: $W_o, W_l: $W_l" )
                        end
                        
                        if pol_c[a_idx, h_idx, z_idx, t]<0
                            c_temp = pol_c[a_idx, h_idx, z_idx, t]
                            println("warning! negative consumption, c = $c_temp at A_i=$a_idx, H_i=$h_idx, Z_i=$z_idx, t=$t")
                        end
                    end
                end
            end
            
            # 다음 연령에서 사용할 interpolation 만들기
            V_itp = Dict()
            Threads.@threads for k=1:K
                Threads.@threads for j=1:J
                interp = LinearInterpolation(a_grid, V[:,j,k,t])
                V_itp["func$(j)_$(k)"] = interp
                end
            end
            
        end
        
        
        # 시간측정 end
    end


    ###################(Step 2 : Stationary distribution) #################
    # transition matrix

    @time begin
        NJK = N*J*K
        W = spzeros(NJK*T, NJK*T)
        stat_mat = repeat(stat_dist', outer=(5, 1))

        # 각 블록을 적절한 위치에 배치
        Threads.@threads for t in 1:T
            rows = ((t-1)*NJK+1):(t*NJK)
            # 죽음
            W[rows, 1:NJK] = BTM(pol_a[:,:,:,t], pol_h[:,:,:,t], stat_mat, N,J,K, a_grid, h_grid).*(1-s_prob[t]/(1+n))  # 그대로 물려받으나 생산성은 랜덤
            # W[rows, 1:NJK] = BTM(zeros(N,J,K), zeros(N,J,K), stat_mat, N,J,K, a_grid, h_grid).*(1-s_prob[t]/(1+n))      # a=0, h=0에서 시작
            # W[rows, 1:NJK] = BTM(ones(N,J,K), zeros(N,J,K), stat_mat, N,J,K, a_grid, h_grid).*(1-s_prob[t]/(1+n))      # a=1, h=0에서 시작 (A는 1/n)
            # 생존
            if t<T
                W[rows, (t*NJK+1):(t*NJK+NJK)] = BTM(pol_a[:,:,:,t], pol_h[:,:,:,t], Pi, N,J,K, a_grid, h_grid).*s_prob[t]/(1+n)   
            end
        end
    
    
        # psi 찾기
        psi_ini = ones(1, NJK*T)
        psi_ini = psi_ini./sum(psi_ini)
        
        max_iter = 5000
        tol2     = 1e-10
        stat_iter = 0
        psi_old = psi_ini
        diff_dist = 10
        
        for i=1:max_iter
            psi_new = psi_old*W
            psi_new = psi_new./sum(psi_new)
            diff_dist = maximum(abs.(psi_new.-psi_old))
            stat_iter +=1
            psi_old = psi_new
            if diff_dist < tol2
                break
            end
        end    

        println("반복횟수: ", stat_iter, "오차: ", diff_dist)
    end
    
        # Aggregation
        psi = psi_old./sum(psi_old)
        psi_mat = reshape(psi,N,J,K,T)
        share_R = sum((pol_h.== 0) .* psi_mat) 
        share_L = sum((pol_h.>pol_s) .* psi_mat)
        share_HRET = sum((q.*pol_h.>cutoff).*psi_mat)
        demand_rental = sum((pol_h.== 0).*pol_s.* psi_mat) 
        supply_rental = sum((pol_h.>pol_s) .*(pol_h-pol_s).* psi_mat)
        excess_demand_rental = demand_rental - supply_rental 
        trading = sum((pol_h.!=H_grid).* psi_mat)
        Retire = zeros(N,J,K,T)
        Retire[:,:,:,R:end] .=1

        # tax
        income_tax_revenue = sum(income_tax.(Y_grid .+ r.*A_grid.*(A_grid.>0)).*psi_mat)
        tau_R_rate = (income_tax.(Y_grid .+ r.*A_grid.*(A_grid.>0))./(Y_grid .+ r.*A_grid.*(A_grid.>0))).*(pol_h.>pol_s)
        rental_tax_revenue = sum(tau_R_rate.*rho.*(pol_h.-pol_s).*(pol_h.>pol_s).*psi_mat)
        property_tax_revenue = sum(tau_p.*q.*pol_h.*psi_mat)
        HRET_tax_revenue = sum((HRET.*(q.*pol_h .- cutoff)).*(q.*pol_h.>=cutoff).*psi_mat)
        pension_expenditure = sum(pen.*Retire.*psi_mat)
        G_BC = (income_tax_revenue + rental_tax_revenue + property_tax_revenue + HRET_tax_revenue) - (Tr_r*share_R + Tr + pension_expenditure)
        
        # earning gini
        Y_mat = zeros(K,T)
        psi_mat2 = zeros(K,T)
        for z_i=1:K
            for t_i=1:T
                Y_mat[z_i,t_i]=mean(Y_grid[:,:,z_i,t_i])
                psi_mat2[z_i,t_i]=sum(psi_mat[:,:,z_i,t_i])
            end
        end
        earning_gini = Gini(vec(psi_mat2),vec(Y_mat),0)

        # asset gini
        net_wealth = zeros(N,J)
        psi_matNW = zeros(N,J)
        for n_i=1:N
            for h_i=1:J
                net_wealth[n_i,h_i] = a_grid[n_i]+q*h_grid[h_i]
                psi_matNW[n_i,h_i] = sum(psi_mat[n_i,h_i,:,:])
                if net_wealth[n_i,h_i]<0 && psi_matNW[n_i,h_i]>0
                    println("warning negative net wealth")
                end
            end
        end
        net_wealth_gini = Gini(vec(psi_matNW),vec(net_wealth),tol2)
        agg_house_demand = sum(pol_h.*psi_mat)

        println("-------------------(parameters )-----------------")
        @printf("rho                           : %.8f\n", rho)
        @printf("delta                         : %.4f\n", delta)
        @printf("phi                           : %.4f\n", phi)
        @printf("cutoff                        : %.4f\n", cutoff)
        @printf("gamma                         : %.4f\n", gamma)
        @printf("b_bar                         : %.4f\n", b_bar)     
        println("-------------------(Targeted Momentum)-----------------")
        @printf("Renter   share                : %.3f\n", share_R)
        @printf("Landlord share                : %.3f\n", share_L)
        @printf("HRET     share                : %.3f\n", share_HRET)
        @printf("trading  share                : %.3f\n", trading)
        @printf("Earning Gini                  : %.3f\n", earning_gini)
        @printf("Wealth Gini                   : %.3f\n", net_wealth_gini)
        println("-------------------(Market clearing)-------------------")
        @printf("House supply                  : %.5f\n", agg_house_demand)
        @printf("government BC error           : %.5f\n", G_BC)
        @printf("rental market clearing error  : %.5f\n", excess_demand_rental)
        println("---------------------------------------------------------")

        target_momentum = [share_R, share_L, share_HRET, trading, net_wealth_gini]

        return target_momentum, excess_demand_rental, G_BC, agg_house_demand, V, pol_a, pol_h, pol_s, pol_c, psi
end

#############################################################
#############################################################
#############################################################

# Load the micro data
Dir = "C:/Users/82106/Documents/4. 2024학년도/HRET/Model/"
Life_mat = readdlm("$(Dir)/Life_mat.txt")
Idio_mat = readdlm("$(Dir)/Idio_mat.txt")

println("------Housing model 추정을 시작합니다.-----------")
start_str = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
println("시작시각: ", start_str)

function calibration(x,grad)
    println("-----------------------")
    println(x)
    q = 1.5
    HRET    = 0.010;       # HRET tax rate (nominal rate)
    h_grid = [0, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 8, 10, 11, 12, 13]
    # h_grid = [0, 2, 3]
    Tr_r  = 0.062      # lump-sum transfer for renters
    Tr    = 0.0279
    gamma = 0.0065
    phi   = 0.0072*q    
    phi_0 = 0.000
    delta = 0.0185
    b_bar   = 0.31
    cutoff  = 12.5*q
  
    # Housing_PE_LC(q, rho, Tr, phi, phi_0, gamma, delta, h_grid, cutoff, b_bar, HRET, Tr_r, ife_mat, Idio_mat)
    target_momentum, ex, g, h_d, V, pol_a, pol_h, pol_s, pol_c, psi = Housing_PE_LC(q, x[1], Tr, x[3]*q, phi_0, x[5], x[2], h_grid, x[4]*q, x[6], HRET, Tr_r, Life_mat, Idio_mat)
    
    # target momenutm
    y = zeros(6)
    y[1] = (ex - 0)*3
    y[2] = (target_momentum[1]    - 0.333)
    y[3] = (target_momentum[2]    - 0.102)
    y[4] = (target_momentum[3]    - 0.010)
    y[5] = target_momentum[4]    - 0.061
    y[6] = target_momentum[5]    - 0.653 # housing 0.707/ wealth 0.653
     # y8 = (g - 0)*2

    writedlm("$(Dir)/value.txt", V, '\t')
    writedlm("$(Dir)/pol_a.txt", pol_a, '\t')
    writedlm("$(Dir)/pol_h.txt", pol_h, '\t')
    writedlm("$(Dir)/pol_s.txt", pol_s, '\t')
    writedlm("$(Dir)/pol_c.txt", pol_c, '\t')
    writedlm("$(Dir)/psi.txt", psi, '\t')
    writedlm("$(Dir)/eq.txt", [q, x, Tr, HRET, gamma, phi, phi_0, delta, b_bar, cutoff, Tr_r], '\t')

    return y
 end
 
 x0 = [0.09970654, 0.019, 0.0076, 12.5, 0.007, 0.30]

 #      rho ,   Tr ,  phi_0,  delta
 # x0 = [0.120, 0.049, 0.00, 0.005]
 
#최적화 문제 설정
opt = Opt(:LN_NELDERMEAD, 6)
#opt = Opt(:LN_COBYLA, 4)
# ftol과 xtol 설정
#lower_bounds!(opt, [0.0900, 0.010, 0.005, 12.9, 0.006, 0.350])
#upper_bounds!(opt, [0.1100, 0.030, 0.010, 12.1, 0.008, 0.250])
opt.initial_step = [0.001, 0.001, 0.001, 0.1, 0.001, 0.01]
opt.ftol_rel = 1e-4
opt.xtol_rel = 1e-4
min_objective!(opt, (x, grad) -> sum(abs2, calibration(x, grad)))

# 최적화 실행
(minf, minx, ret) = NLopt.optimize(opt, x0)

println("해: ", minx)
println("오차: ", minf)

# calibration(0.352)
# sol = find_zero(calibration, 0.09826760)
# sol = find_zero(calibration, (0.100, 0.250), Roots.Brent())
# println("해: ", sol)
# result = nlsolve(calibration, x0, ftol = 1e-4, xtol = 1e-4)
# 결과 출력
# println("해: ", result.zero)
#println("수렴 여부: ", result.converged)

end_str = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
println("종료시각: ", end_str)