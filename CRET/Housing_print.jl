# Housing model with a Life cycle 
# Yeongwoong Do (2024-04-26)

using LinearAlgebra, Interpolations, Plots, Printf, Dates, MAT, Optim, SparseArrays, NLopt, Roots, DelimitedFiles, Statistics

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
function housing_tax_cost(gamma, delta, h, h_prime, h_adj, cutoff, q, tau_p, HRET)
    transaction_cost = gamma*(q*h_prime + q*(1-delta)*h)*(h_prime!=h)
    property_tax = tau_p*q*h_prime
    HRET_tax = HRET*h_adj*(q*h_prime - cutoff)*(q*h_prime >= cutoff)
    return transaction_cost + property_tax + HRET_tax
end

# RHS of Renter's problem
function RHS_renter(a_prime, bar_x, sigma, alpha, beta, V_itp, Pi, s_prob, rho, z_idx, t, b_bar, gamma, q, h, delta) 
    x = bar_x - a_prime -gamma*q*(1-delta)*h
    if x <0
        return -1e6
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
        return -1e6
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
function RHS_landlord(a_prime, h_prime, h_prime_idx, bar_x, sigma, alpha, beta, V_itp, Pi, s_prob, rho, q, z_idx, t, r_adj, tau_R, h_cost, b_bar, phi)
    x = bar_x -a_prime - h_cost - (q - rho*(1-tau_R*r_adj)+ phi)*h_prime
    s = (1-alpha)*x/(rho*(1-tau_R*r_adj)-phi)
    if x <0         # c> 0, s>0
        return -1e6
    end
    if h_prime <= s  # h'>s
        return -1e6
    end
    T = length(s_prob)
    u = (x^(1-sigma)/(1-sigma))*(1/(rho*(1-tau_R*r_adj)))^((1-alpha)*(1-sigma))
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
function Gini(pop,x)
    # pop의 원소가 0보다 큰 인덱스만 필터링
    # 이것을 하지 않으면 아무도 비싼 자산을 안가지고 있는데, 가지고 있는 것처럼 나와서 지니계수가 이상하게 나옴
    valid_idx = findall(p -> p > 1e-10, pop)
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

function Housing_PE_LC(q, h_min, rho, phi, cutoff, gamma, b_bar, delta, J, HRET, Life_mat, Idio_mat)
    @printf("rho: %.8f\n", rho)
    # save vector
    w_age  = Life_mat[:,1]
    T = length(w_age)
    s_prob = Life_mat[:,2]
    mu     = Life_mat[:,3]
    z_grid = Idio_mat[:,1]
    K = length(z_grid)
    Pi     = Idio_mat[:,2:end]

    # parameters
    sigma   = 2.000;       # inverse of EIS
    alpha   = 0.778;       # share of non-housing consumption
    r       = 0.024;       # interest rate
    kappa   = 0.013;       # spread of mortgage debt
    beta    = 0.975;       # time discount factor
    LTV     = 0.600;       # LTV limit
    tau_p   = 0.001;       # property tax rate (effective rate)
    tau_R   = 0.140;       # rental income tax rate (effective rate)
    r_adj   = 0.500;       # adjustment coefficient for rental income tax
    h_adj   = 0.600;        # adjustment coefficient for HRET
    pen     = 0.149;       # pension 
    R       = 61-19;       # retire age
    DTI     = 0.500;        # DTI limit
    n       = 0;           # net population growth rate 
    minus_inf = -1e6;

    # A grid, H_grid
    a_min = -10
    a_max = 20
    N = 100
    h_max = 15
    # a_grid 생성
    a_debt_grid = -(LinRange(sqrt(-a_min), sqrt(0), Int(N/2))).^2
    delta = a_debt_grid[end] - a_debt_grid[end-1]
    a_saving_grid = exp.(LinRange(log(delta), log(a_max), Int(N/2)))
    a_grid = [a_debt_grid; a_saving_grid]
    a_lb_renter = 0
    # h_grid 생성
    h_grid = range(h_min, h_max, J-1)
    h_grid = [0; h_grid];
    println(h_grid)

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
    
    for t=T:-1:1  # Backward!  
            
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

                        bar_x = y + (1+r)*a  + kappa*a*(a<0) + q*(1-delta)*h -income_tax(tax_base) + pen*(t>=R)    # maximum expenditure

                        ###########################################
                        # Renter's problem W(.,h'=0)
                        
                        a_ub_renter = min(bar_x-gamma*q*(1-delta)*h , a_max)
                        if a_ub_renter < a_lb_renter # 말도 안되는 상황이 오면
                            W_R = minus_inf
                            W_R_pol_a = 0
                            W_R_pol_s = 0 
                            W_R_pol_c = 0 
                        else
                            function f_rent(x)
                                return -RHS_renter(x, bar_x, sigma, alpha, beta, V_itp, Pi, s_prob, rho, z_idx, t, b_bar, gamma, q, h, delta) 
                            end
                            results = Optim.optimize(f_rent, a_lb_renter, a_ub_renter)
                            W_R = -Optim.minimum(results)
                            W_R_pol_a = Optim.minimizer(results)
                            W_R_pol_s = (1-alpha)*(bar_x- W_R_pol_a)/rho
                            W_R_pol_c = alpha*(bar_x- W_R_pol_a)
                        end
                        ################################################
                        
                        ############################################
                        # Owner-occupier's problem W(.,h'=s)
                        W_o_temp = zeros(J-1);
                        pol_o_a_temp = zeros(J-1);
            
                        Threads.@threads for h_prime_idx = 2:J
                            h_prime = h_grid[h_prime_idx]
                            h_cost = housing_tax_cost(gamma, delta, h, h_prime, h_adj, cutoff, q, tau_p, HRET)
                            if t == T
                                a_lb_owner = 0  # t==T일 경우, 0을 사용
                            else
                                a_lb_owner = max(-LTV*q*h_prime, -DTI*y/(0.1+r+kappa), a_min)
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
                        W_o_pol_c = bar_x - W_o_pol_a - q*W_o_pol_h - housing_tax_cost(gamma, delta, h, W_o_pol_h, h_adj, cutoff, q, tau_p, HRET)

                        ################################################################

                        ################################################################
                        # Landlord's problem W(.,h'>s)
                        W_l_temp = zeros(J-1);
                        pol_l_a_temp = zeros(J-1);
                        pol_l_s_temp = zeros(J-1);   

                        Threads.@threads for h_prime_idx = 2:J
                            h_prime = h_grid[h_prime_idx]
                            h_cost = housing_tax_cost(gamma, delta, h, h_prime, h_adj, cutoff, q, tau_p, HRET)
                            if t == T
                                a_lb_landlord = 0  # t==T일 경우, 0을 사용
                            else
                                a_lb_landlord = max(-LTV*q*h_prime, -DTI*y/(0.1+r+kappa), a_min)
                            end
                            a_ub_landlord = min(bar_x - q*h_prime - h_cost, a_max)

                            if a_ub_landlord < a_lb_landlord # 말도 안되는 상황이 오면
                                W_l_temp[h_prime_idx-1]=minus_inf
                                pol_l_a_temp[h_prime_idx-1]= 0
                                pol_l_s_temp[h_prime_idx-1]= 0
                            else
                                function f_landlord(x)
                                    return - RHS_landlord(x, h_prime, h_prime_idx, bar_x, sigma, alpha, beta, V_itp, Pi, s_prob, rho, q, z_idx, t, r_adj, tau_R, h_cost, b_bar, phi)
                                end
                                results_landlord = Optim.optimize(f_landlord, a_lb_landlord, a_ub_landlord)
                                W_l_temp[h_prime_idx-1] = -Optim.minimum(results_landlord)
                                pol_l_a_temp[h_prime_idx-1] = Optim.minimizer(results_landlord)
                                pol_l_s_temp[h_prime_idx-1] = (1-alpha)*(bar_x -pol_l_a_temp[h_prime_idx-1] - h_cost - (q-rho*(1-tau_R*r_adj)+phi)*h_prime)/(rho*(1-tau_R*r_adj)-phi)
                            end

                        end

                        W_l = maximum(W_l_temp)
                        l_idx = argmax(W_l_temp)
                        W_l_pol_a = pol_l_a_temp[l_idx]
                        W_l_pol_h = h_grid[l_idx+1]
                        W_l_pol_s = pol_l_s_temp[l_idx]
                        W_l_pol_c = W_l_pol_s*(rho*(1-tau_R*r_adj))*alpha/(1-alpha)

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
    
        # 각 블록을 적절한 위치에 배치
        Threads.@threads for t in 1:T
            rows = ((t-1)*NJK+1):(t*NJK)
            # 죽음
            # W[rows, 1:NJK] = BTM(pol_a[:,:,:,t], pol_h[:,:,:,t], Pi, N,J,K, a_grid, h_grid).*(1-s_prob[t]/(1+n))  # 그대로 물려받음
            W[rows, 1:NJK] = BTM(zeros(N,J,K), zeros(N,J,K), Pi, N,J,K, a_grid, h_grid).*(1-s_prob[t]/(1+n))      # a=0, h=0에서 시작
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
        share_R = sum((pol_stance.== 0) .* psi_mat) 
        share_L = sum((pol_stance.== 2) .* psi_mat)
        share_HRET = sum((q.*pol_h.>cutoff).*psi_mat)
        demand_rental = sum((pol_stance.== 0).*pol_s.* psi_mat) 
        supply_rental = sum((pol_stance.== 2).*(pol_h-pol_s).* psi_mat)
        excess_demand_rental = demand_rental - supply_rental 
        
        # inequality
        psi = vec(psi)
        total_income = Y_grid .+ r.*A_grid.*(A_grid.>0) .+ rho.*(pol_h.-pol_s).*(pol_h.>pol_s)
        Retire = zeros(N,J,K,T)
        Retire[:,:,:,R:end] .=1
        after_tax_income = total_income.-income_tax.(Y_grid .+ r.*A_grid.*(A_grid.>0)) .-tau_R.*r_adj.*rho.*(pol_h.-pol_s).*(pol_h.>pol_s) .+ pen.*Retire .-tau_p.*q.*pol_h .- HRET.*h_adj.*(q.*pol_h .- cutoff).*(q.*pol_h .>= cutoff)
        earning_gini = Gini(psi,vec(Y_grid)) 
        income_gini = Gini(psi,vec(total_income))
        desposible_gini = Gini(psi,vec(after_tax_income))
        net_wealth_gini = Gini(psi,vec(A_grid .+ q.*H_grid ))
        housing_asset_gini = Gini(psi,vec(H_grid))
        consumption_gini = Gini(psi,vec(pol_c))

        # housing and mortgage
        trading = sum((pol_h.!=H_grid).* psi_mat)
        homeowner_end = sum((pol_h[:,:,:,end-9:end-5].> 0).* psi_mat[:,:,:,end-9:end-5])/sum(psi_mat[:,:,:,end-9:end-5])
        indebted = sum((pol_a.<0).* psi_mat)
        LTV = (-pol_a)./(q.*pol_h)
        avg_LTV = mean(LTV[pol_a.<0])

        # tax
        income_tax_revenue = sum(income_tax.(Y_grid .+ r.*A_grid.*(A_grid.>0)).*psi_mat)
        rental_tax_revenue = sum(tau_R.*r_adj.*rho.*(pol_h.-pol_s).*(pol_h.>pol_s).*psi_mat)
        property_tax_revenue = sum(tau_p.*q.*pol_h.*psi_mat)
        HRET_tax_revenue = sum(HRET.*h_adj.*(q.*pol_h .- cutoff).*(q.*pol_h .>= cutoff).*psi_mat)
        pension_expenditure = sum(pen.*Retire.*psi_mat)
        exogenous_expenditure = income_tax_revenue + rental_tax_revenue + property_tax_revenue + HRET_tax_revenue - pension_expenditure

        # HRET
        share_HRET_landlord = sum((q.*pol_h.>cutoff).*(pol_h.>pol_s).*psi_mat)
        share_HRET_61_plus = sum((q.*pol_h.>cutoff).*Retire.*psi_mat)
        HRET_property = HRET_tax_revenue/property_tax_revenue


        println("-------------------(Targeted Momentum)-----------------")
        @printf("Renter   share                : %.3f\n", share_R)
        @printf("Landlord share                : %.3f\n", share_L)
        @printf("HRET     share                : %.3f\n", share_HRET)
        @printf("Earning Gini                  : %.3f\n", earning_gini)
        println("------------(Untargeted Momentum; inequality)-----------")
        @printf("Income Gini                   : %.3f\n", income_gini)
        @printf("desposible Gini               : %.3f\n", desposible_gini)
        @printf("Net wealth Gini               : %.3f\n", net_wealth_gini)
        @printf("Housing Gini                  : %.3f\n", housing_asset_gini)
        @printf("consumption Gini              : %.3f\n", consumption_gini)
        println("------------(Untargeted Momentum; housing )-----------")
        @printf("House supply                  : %.3f\n", sum(pol_h.*psi_mat))
        @printf("Aggregate shelter             : %.3f\n", sum(pol_s.*psi_mat))
        @printf("trading  share                : %.3f\n", trading)
        @printf("indebted household            : %.3f\n", indebted)
        @printf("average LTV                   : %.3f\n", avg_LTV)
        @printf("75~80 homeowner share         : %.3f\n", homeowner_end)
        @printf("excess demand in rental market: %.3f\n", excess_demand_rental)
        println("-------------------(Government )-----------------------")
        @printf("Income tax revenue            : %.3f\n", income_tax_revenue)
        @printf("Rental income tax revenue     : %.3f\n", rental_tax_revenue)
        @printf("Property tax revenue          : %.3f\n", property_tax_revenue)
        @printf("HRET tax revenue              : %.3f\n", HRET_tax_revenue)
        @printf("Pension expnediture           : %.3f\n", pension_expenditure)
        @printf("Exogenous expenditure(g)      : %.3f\n", exogenous_expenditure)
        println("-------------------(HRET )-----------------------")
        @printf("Landlord fraction             : %.3f\n", share_HRET_landlord/share_HRET)
        @printf("61 + fraction                 : %.3f\n", share_HRET_61_plus/share_HRET)
        @printf("HRET / Property tax           : %.3f\n", HRET_property)
        println("---------------------------------------------------------")

        target_momentum = [share_R, share_L, share_HRET, trading, homeowner_end]

        return target_momentum, excess_demand_rental, V, pol_a, pol_h, pol_s, pol_c, psi
end

#############################################################
#############################################################
#############################################################

# Load the micro data
file = matopen("C:/Users/82106/Documents/4. 2024학년도/HRET/calibration/Life_mat.mat")
Life_mat = read(file, "Life_mat")
close(file)
file = matopen("C:/Users/82106/Documents/4. 2024학년도/HRET/calibration/Idio_mat.mat")
Idio_mat = read(file, "Idio_mat")
close(file)

println("------Housing model 추정을 시작합니다.-----------")
start_str = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
println("시작시각: ", start_str)

# HRET =0.01  q= 2, rho = 0.10994272648006477

HRET = 0.0
q = 2.00384986
h_min = 1.415
rho = 0.10993888
phi = 0.01500
cutoff = 21
gamma = 0.02
b_bar = 0.50
delta = 0.03
J     = 30

println(h_min, rho, phi, delta)

target_momentum, ex, V, pol_a, pol_h, pol_s, pol_c, psi = Housing_PE_LC(q, h_min, rho, phi, cutoff, gamma, b_bar, delta, J, HRET, Life_mat, Idio_mat)

filepath1 = raw"C:\\Users\\82106\\Documents\\4. 2024학년도\\HRET\\Model\\value.txt"
filepath2 = raw"C:\\Users\\82106\\Documents\\4. 2024학년도\\HRET\\Model\\pol_a.txt"
filepath3 = raw"C:\\Users\\82106\\Documents\\4. 2024학년도\\HRET\\Model\\pol_h.txt"
filepath4 = raw"C:\\Users\\82106\\Documents\\4. 2024학년도\\HRET\\Model\\pol_s.txt"
filepath5 = raw"C:\\Users\\82106\\Documents\\4. 2024학년도\\HRET\\Model\\pol_c.txt"
filepath6 = raw"C:\\Users\\82106\\Documents\\4. 2024학년도\\HRET\\Model\\psi.txt"

writedlm(filepath1, V, '\t')
writedlm(filepath2, pol_a, '\t')
writedlm(filepath3, pol_h, '\t')
writedlm(filepath4, pol_s, '\t')
writedlm(filepath5, pol_c, '\t')
writedlm(filepath6, psi, '\t')

end_str = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
println("종료시각: ", end_str)