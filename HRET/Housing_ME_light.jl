# Housing model with a Life cycle 
# Yeongwoong Do (2024-04-26)

using LinearAlgebra, Interpolations, Plots, Printf, Dates, MAT, Optim, SparseArrays, NLopt, Roots, DelimitedFiles, Statistics, NLsolve, ProgressMeter

Dir = "C:/Users/82106/Documents/4. 2024학년도/HRET/Model/"
include("$(Dir)/help_function.jl")
using .help_function

include("$(Dir)/data_params_grid.jl")
using .data_params_grid

#############################################################
####################(Main)###################################
#############################################################

function Housing_PE_LC(q, rho, Tr, HRET, tau_tr, earning_info, grid)
    
    # earning information
    w_age       = earning_info[1]
    s_prob      = earning_info[2]
    mu          = earning_info[3]
    z_grid      = earning_info[4]
    Pi          = earning_info[5]
    stat_dist   = earning_info[6]
    T           = earning_info[7]
    K           = earning_info[8] 

    # calibrated parameters
    gamma = 0.007596959553673047
    phi   = 0.007672072358547014*1.5    
    delta = 0.01912987211359198
    b_bar   = 0.3027267224200313 
    cutoff  = 12.498110381496865*1.5
    H_s     = 2.943466566628909
    
    # parameters
    sigma   = 2.000;       # inverse of EIS
    alpha   = 0.778;       # share of non-housing consumption
    r       = 0.024;       # interest rate
    kappa   = 0.013;       # spread of mortgage debt
    beta    = 0.975;       # time discount factor
    LTV     = 0.600;       # LTV limit
    tau_p   = 0.001;       # property tax rate (effective rate)
    pen     = 0.149;       # pension 
    phi_0   = 0.000        # landlord fixed cost
    Tr_r    = 0.062      # lump-sum transfer for renters
    R       = 37;          # retire age
    n       = 0;
    minus_inf = -1e10;

    # A grid, H_grid
    a_grid  = grid[1]
    h_grid  = grid[2]
    N       = grid[3]
    J       = grid[4]
    
    # fast test grid
    # a_grid = [-1,0,1]
    # h_grid = [0,1,2]
    # N = length(a_grid)
    # J = length(h_grid)

    a_lb_renter = 0
    a_min = a_grid[1]
    a_max = a_grid[N]

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
                            h_cost = housing_tax_cost(gamma, delta, h, h_prime, cutoff, q, tau_p, HRET, tau_tr)
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
                        W_o_pol_c = bar_x - W_o_pol_a - q*W_o_pol_h - housing_tax_cost(gamma, delta, h, W_o_pol_h, cutoff, q, tau_p, HRET, tau_tr)

                        ################################################################

                        ################################################################
                        # Landlord's problem W(.,h'>s)
                        W_l_temp = zeros(J-1);
                        pol_l_a_temp = zeros(J-1);
                        pol_l_s_temp = zeros(J-1);   

                        Threads.@threads for h_prime_idx = 2:J
                            h_prime = h_grid[h_prime_idx]
                            h_cost = housing_tax_cost(gamma, delta, h, h_prime, cutoff, q, tau_p, HRET, tau_tr)
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
                            println("error ! landlord contradiction, H_prime :$W_l_pol_h, S: $W_l_pol_s, W_l: $W_l")
                        end

                        ###################################        
                        
                        max_temp = max(W_R, W_o, W_l)

                        if W_R == max_temp # HH will be a renter
                            V[a_idx, h_idx, z_idx, t]=W_R
                            pol_stance[a_idx, h_idx, z_idx, t]=0
                            pol_a[a_idx, h_idx, z_idx, t] = W_R_pol_a
                            pol_h[a_idx, h_idx, z_idx, t] = 0
                            pol_s[a_idx, h_idx, z_idx, t] = W_R_pol_s
                            pol_c[a_idx, h_idx, z_idx, t] = W_R_pol_c
                        elseif W_o == max_temp # HH will be a owner-occupier
                            V[a_idx, h_idx, z_idx, t]=W_o
                            pol_stance[a_idx, h_idx, z_idx, t]=1
                            pol_a[a_idx, h_idx, z_idx, t] = W_o_pol_a
                            pol_h[a_idx, h_idx, z_idx, t] = W_o_pol_h
                            pol_s[a_idx, h_idx, z_idx, t] = W_o_pol_s
                            pol_c[a_idx, h_idx, z_idx, t] = W_o_pol_c
                        elseif W_l == max_temp # HH will be a landlord
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
        demand_rental = sum((pol_h.== 0).*pol_s.* psi_mat) 
        supply_rental = sum((pol_h.>pol_s) .*(pol_h-pol_s).* psi_mat)
        excess_demand_rental = demand_rental - supply_rental 
        H_d = sum(pol_h.*psi_mat)
        excess_demand_housing = H_d - H_s

        # Tax
        Retire = zeros(N,J,K,T)
        Retire[:,:,:,R:end] .=1
        share_R = sum((pol_h.== 0).* psi_mat) 
        income_tax_revenue = sum(income_tax.(Y_grid .+ r.*A_grid.*(A_grid.>0)).*psi_mat)
        tau_R_rate = (income_tax.(Y_grid .+ r.*A_grid.*(A_grid.>0))./(Y_grid .+ r.*A_grid.*(A_grid.>0))).*(pol_h.>pol_s)
        rental_tax_revenue = sum(tau_R_rate.*rho.*(pol_h.-pol_s).*(pol_h.>pol_s).*psi_mat)
        property_tax_revenue = sum(tau_p.*q.*pol_h.*psi_mat)
        HRET_tax_revenue = sum((HRET.*(q.*pol_h .- cutoff)).*(q.*pol_h.>=cutoff).*psi_mat)
        pension_expenditure = sum(pen.*Retire.*psi_mat)
        mansion_tax_revenue = sum(tau_tr.*(q.*pol_h.*(q.*pol_h.>=cutoff) .+ q.*(1-delta).*H_grid.*(q.*(1-delta).*H_grid.>=cutoff)).*(pol_h.!=H_grid).*psi_mat)
        total_tax_revenue = income_tax_revenue + rental_tax_revenue + property_tax_revenue + HRET_tax_revenue + mansion_tax_revenue
        G_BC = total_tax_revenue - (Tr_r*share_R + Tr + pension_expenditure)

        println("----------------(Price and Transfer)-----------------")
        @printf("q                             : %.8f\n", q)
        @printf("rho                           : %.8f\n", rho)
        @printf("Tr                            : %.8f\n", Tr)
        @printf("HRET                          : %.4f\n", HRET)
        println("-------------------(Market clearing)-------------------")
        @printf("Housing marekt clearing error : %.5f\n", excess_demand_housing)
        @printf("rental market clearing error  : %.5f\n", excess_demand_rental)
        @printf("government BC error           : %.5f\n", G_BC)
        println("---------------------------------------------------------")

        return excess_demand_rental, G_BC, excess_demand_housing, V, pol_a, pol_h, pol_s, pol_c, psi
end

#############################################################
#############################################################
#############################################################

# Load the micro data
earning  = micro_data("$(Dir)")
grid   = a_h_grid("$(Dir)")

println("------Housing model 추정을 시작합니다.-----------")
start_str = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
println("시작시각: ", start_str)

function obj(x)
    HRET    = 0.020      # HRET tax rate (nominal rate)
    tau_tr  = 0.000       # mansion tax rate
    ex_r, g, ex_h, V, a, h, s, c, psi = Housing_PE_LC(x[1], x[2], x[3], HRET, tau_tr, earning, grid)
    # market clearing and GBC condition
    y = zeros(3)
    y[1] = (ex_r - 0)
    y[2] = (g - 0)
    y[3] = (ex_h - 0)

    writedlm("$(Dir)/value.txt", V, '\t')
    writedlm("$(Dir)/pol_a.txt", a, '\t')
    writedlm("$(Dir)/pol_h.txt", h, '\t')
    writedlm("$(Dir)/pol_s.txt", s, '\t')
    writedlm("$(Dir)/pol_c.txt", c, '\t')
    writedlm("$(Dir)/psi.txt", psi, '\t')
    writedlm("$(Dir)/eq.txt", [x, HRET], '\t')

    return y
 end
 
#       q ,      rho ,              Tr
x0 = [1.499150, 0.0996850, 0.027710]

result = nlsolve(obj, x0, ftol = 1e-4, xtol = 1e-4)
# 결과 출력
println("해: ", result.zero)
println("수렴 여부: ", result.converged)

end_str = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
println("종료시각: ", end_str)