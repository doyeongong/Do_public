# New version of HRET(High-value Real Estate Tax)
# edited by Yeongwoong Do (2024.09.05)
# Transition path

using LinearAlgebra, Interpolations, Plots, Printf, Dates, MAT
using Optim, SparseArrays, NLopt, Roots, DelimitedFiles, Statistics, NLsolve, ProgressMeter


Dir = "C:/Users/82106/Documents/4. 2024학년도/HRET/Model"

include("$(Dir)/A1_Rouwenhorst.jl")
include("$(Dir)/A2_Grid.jl")
include("$(Dir)/A3_transition_matrix.jl")
include("$(Dir)/A4_root_finding.jl")
include("$(Dir)/A5_Korean_tax.jl")
include("$(Dir)/A6_Inequality.jl")

include("$(Dir)/B1_Household.jl")
include("$(Dir)/B2_Fitting.jl")

using .Rouwenhorst
using .Grid
using .transition_matrix
using .root_finding
using .Korean_tax
using .Inequality

using .Household
using .fitting

##############  1. Life cycle earnings ##############
# Load the micro data
Life_mat = readdlm("$(Dir)/Life_mat.txt")
Idio_mat = readdlm("$(Dir)/Idio_z3.txt",',')

# Life cycle (SFLC and Korean statistics)
w_age        = Life_mat[:,1]
s_prob       = Life_mat[:,2]
s_prob[end]  = 0
mu           = Life_mat[:,3]
w_age        = w_age./(w_age'*mu)
T            = length(w_age)
K, _         = size(Idio_mat)
z_grid       = Idio_mat[:,1]
Π            = Idio_mat[:,2:K+1]
stat_dist    = Idio_mat[:,end];

##############  2. parameters #######################
# parameters

param = Dict(
        # Life cycle
        "T"       =>  T,              # Total period of life time 
        "R"       =>  37,             # retirement (starting point of get a pension)    
        "p"       =>  0,              # pension / mean earnings
        # Grid
        "N_z"     =>  K,              # Points in z grid
        "N_a"     =>  120,            # Points in mu grid
        "a_min"   =>  -20,            # lower bound for mu
        "a_max"   =>  30,             # upper bound for mu
        "N_J"     =>  8,              # Points in h grid
        "h_min"   =>  1.22,           # lower bound for h
        "h_max"   =>  10.0,           # upper bound for h
        # External calibration
        "β"       =>  0.970,          # time discount factor
        "r"       =>  0.024,          # deposit interest rate (BOK)
        "κ"       =>  0.013,          # spread of mortgage debt (BOK)
        "θ"       =>  0.700,          # LTV limit
        "DTI"     =>  0.600,          # DTI limit
        "m"       =>  0.00422,        # moving cost (NaSTaB)
        "δ_R"     =>  0.015,          # share of rental unit cost
        "τ_R"     =>  0.13,           # separately rental income tax
        "unit"    =>  5621,           # mean earnig
        "s_min"   =>  0.00,           # minimum unit of rental housing supply 
        "-∞"      =>  -1e10,              
        # Preference - endogeneous calibration
        "σ"       =>  4.50,          # curvature of C
        "ξ"       =>  1.152,          # curvature of S
        "α"       =>  0.18,          # weight for shelter service in utility
        "b_bar"   =>  4.20,          # weight for bequest in utility
        "ψ"       =>  1.20,          # curvature of bequest
        # Housing - endogeneous calibration
        "ϕ"       =>  0.0150,          #################  landlord fixed cost (landlord target)
        "γ_s"     =>  0.0250,          #################  housing selling cost (trading share target)
        "δ"       =>  0.0317,          # depreciation for housing stock (Cho, Lee and Do(2012, QNAR)) 0.0317
        "adj"     =>  1.100,            #################  housing value adjustment for property tax (HRET target) (낮출수록 내는 비율 낮아짐)
        "H_s"     =>  1.5001241        # H_s
);

##############  3. Grid #####################################################
a_grid         = densezero(param["a_min"], param["a_max"], param["N_a"])
h_grid = [0.0, 1.22, 2.0, 3.0, 4.07, 5.825, 8.0, 10.0] 
# h_grid         = [0.0; range(param["h_min"], param["h_max"], param["N_J"]-1)]
# h_grid         = [0.0; log_type_grid(param["h_min"], param["h_max"], param["N_J"]-1)]
# println(h_grid)
H_grid         = zeros(param["N_a"], param["N_J"], K, T)
for i=1:param["N_J"]
    H_grid[:,i,:,:].=h_grid[i]
end
Y_grid         = zeros(param["N_a"], param["N_J"], K, T)
for k=1:K
    for t=1:T
        Y_grid[:,:,k,t].=w_age[t]*z_grid[k]
    end
end
Grid_param    = param["N_a"], param["N_J"], param["N_z"], param["T"], param["a_min"], param["a_max"], param["unit"], Π, s_prob, h_grid, a_grid;
N, J, K, T, a_min, a_max, unit, Π, s_prob, h_grid, a_grid = Grid_param;
housing_param = param["γ_s"], param["m"], param["r"], param["θ"], param["κ"], param["δ"], param["r"]+param["κ"] , param["δ_R"], param["ϕ"]
γ_s, m, r, θ, κ, δ, r_m, δ_R, ϕ = housing_param
p, R, adj, minus_infi = param["p"], param["R"], param["adj"], param["-∞"]


############## 4. Household problem ############################################

function Housing_PE_LC(q, ρ, Tr, Grid_param, V_next)   
       
    # parameters
    N, J, K, T, a_min, a_max, unit, Π, s_prob, h_grid, a_grid = Grid_param
    util_param    = param["β"], param["α"], param["σ"], param["ξ"], param["b_bar"] , param["ψ"]
    housing_param = param["γ_s"], param["m"], param["r"], param["θ"], param["κ"], param["δ"], param["r"]+param["κ"] , param["δ_R"], param["ϕ"]
    p, R, adj, minus_infi = param["p"], param["R"], param["adj"], param["-∞"]
    γ_s, m, r, θ, κ, δ, r_m, δ_R, ϕ = housing_param
    
    # empty-set
    V, pol_a, pol_c, pol_s, pol_h, total_tax = zeros(N,J,K,T), zeros(N,J,K,T), zeros(N,J,K,T), zeros(N,J,K,T), zeros(N,J,K,T), zeros(N,J,K,T)
    V_itp = Dict()

@time begin
            
    @showprogress for t=T:-1:1  # Backward!
            # age-specific housing tax vector
            # h_tax = (property_tax.(q*h_grid*unit*adj).+ HRET.(q*h_grid*unit*adj,t).*2)*q.*h_grid
            h_tax = (property_tax.(q*h_grid*unit*adj))*q.*h_grid
            # println("h_tax: $(h_tax)")
            
            for z_idx=1:K
                z = z_grid[z_idx]  
                y = w_age[t]*z
                DTI_lim = - param["DTI"]*y/r_m
                
                for h_idx=1:J
                    h = h_grid[h_idx]
                    j_p_lower = Int(2)
                                       
                    for a_idx=1:N
                        a = a_grid[a_idx]        
                        tax_base = y + r*a*(a>=0)
                        τ_y     = income_tax(tax_base*unit)
                        R_de     = rent_deduct(tax_base*unit)
                        common_x = y + (1+r)*a + (1-δ)*q*h + κ*a*(a<0) - τ_y*tax_base + p*(t>=R) + Tr                      
                        
                        W_temp, a_temp, h_temp, s_temp, c_temp = zeros(3), zeros(3), zeros(3), zeros(3), zeros(3)
                        τ_R_in = minimum([τ_y, param["τ_R"]])
                        # τ_R_in = param["τ_R"]
                        
                        ####### (Renter) ##############
                        W_temp[1], a_temp[1], c_temp[1], s_temp[1] = Renter_problem(q, ρ, h, t, z_idx, common_x, 
                                                                V_itp, R_de, util_param, Grid_param, housing_param)
                        
                        if s_temp[1] > param["h_max"]
                            W_temp[1] = minus_infi
                        end
                        
                        ######## (Owner-occupier) #####
                        W_temp[2], a_temp[2], c_temp[2], h_temp[2] = Owner_problem(j_p_lower, q, h, t, z_idx, common_x, 
                                                            V_itp, h_tax, DTI_lim, util_param, Grid_param, housing_param)
                        s_temp[2] = h_temp[2]

                        if W_temp[2] <= W_temp[1] 
                            W_temp[3] = minus_infi
                        else
                            ######## (Landlord) ###########
                            j_p_lower_landlord = findfirst(x -> x == h_temp[2], h_grid)
                            W_temp[3], a_temp[3], c_temp[3], s_temp[3], h_temp[3] = Landlord_problem(j_p_lower, q, ρ, h, t, z_idx, common_x,
                                                            V_itp, h_tax, τ_R_in, DTI_lim, util_param, Grid_param, housing_param)
                            if h_temp[3] - s_temp[3] < param["s_min"]
                                W_temp[3] = minus_infi
                            end
                        end
                        ###############################
                        V[a_idx, h_idx, z_idx, t]       = maximum(W_temp)
                        temp_idx                        = argmax(W_temp)
                        pol_a[a_idx, h_idx, z_idx, t]   = a_temp[temp_idx]
                        pol_c[a_idx, h_idx, z_idx, t]   = c_temp[temp_idx]
                        pol_s[a_idx, h_idx, z_idx, t]   = s_temp[temp_idx]
                        pol_h[a_idx, h_idx, z_idx, t]   = h_temp[temp_idx]                      
                        total_tax[a_idx, h_idx, z_idx, t] = τ_y*tax_base 
                        + τ_R_in*(ρ-δ_R*q)*(h_temp[temp_idx]-s_temp[temp_idx])*(temp_idx==3)
                        # + (property_tax(q*h_temp[temp_idx]*unit*adj) + HRET(q*h_temp[temp_idx]*unit*adj,t).*2)*q*h_temp[temp_idx]
                        + (property_tax(q*h_temp[temp_idx]*unit*adj))*q*h_temp[temp_idx]
                        - R_de*ρ*s_temp[temp_idx]*(temp_idx==1)
                        - mort_credit(a_temp[temp_idx], q, h_temp[temp_idx], unit, r_m)*(temp_idx==2)
                    end
                end
            end
            
            # 다음 연령에서 사용할 interpolation 만들기
            V_itp = Dict()
            for k=1:K
                for j=1:J
                interp = LinearInterpolation(a_grid, V_next[:,j,k,t])
                V_itp["func$(j)_$(k)"] = interp
                end
            end
            
        end
        
        # 시간측정 end
    end
    
    return V, pol_a, pol_h, pol_s, pol_c, total_tax
end

#################### 6. Transition path ###############################

function log_type_guess(seq_start, seq_end, Tran_p)
    distance = abs(seq_start - seq_end)
    if seq_start < seq_end # increasing
        log_grid = -1.0 ./ range(1, stop=Tran_p, length=Tran_p)
        log_grid = log_grid.*distance .+ seq_end
    else # decreasing
        log_grid = 1.0 ./ range(1, stop=Tran_p, length=Tran_p)
        log_grid = log_grid.*distance .+ seq_end
    end
    return log_grid
end

function Transition_one(q_seq, ρ_seq, V_end, D_ini, Tran_p, Grid_param)
    N, J, K, T, a_min, a_max, unit, Π, s_prob, h_grid, a_grid = Grid_param
    NJKT = N*J*K*T
    # reshape V matrix
    V_end = reshape(V_end, N, J, K, T)

    # saving the results
    W_seq = [spzeros(NJKT, NJKT) for _ in 1:Tran_p]
    pol_a_seq = spzeros(NJKT, Tran_p)
    pol_h_seq = spzeros(NJKT, Tran_p)
    pol_s_seq = spzeros(NJKT, Tran_p)
    pol_c_seq = spzeros(NJKT, Tran_p)
    D_seq = zeros(Tran_p+1, NJKT)
    D_seq[1,:] = D_ini
    Hmcc_error = zeros(Tran_p)
    Rmcc_error = zeros(Tran_p)
    gini_seq = zeros(Tran_p,4)
    
    # Backward iteration
    
    V_next = copy(V_end)
    V_tr = copy(V_end)
    for period = Tran_p:-1:1    
        println("Backward iteration: $(period)/$(Tran_p)")
        V, pol_a, pol_h, pol_s, pol_c, _ = Housing_PE_LC(q_seq[period], ρ_seq[period], Tr, Grid_param, V_next)   
        W = tran_Lambda(pol_a, pol_h, Grid_param, stat_dist)
        W_seq[period] = W
        pol_a_seq[:,period] = vec(pol_a)
        pol_h_seq[:,period] = vec(pol_h)
        pol_s_seq[:,period] = vec(pol_s)
        pol_c_seq[:,period] = vec(pol_c)

        # backward update
        V_next = copy(V)
        if period ==1
            V_tr = copy(V)
        end
    end
    println("Backward iteration completed")

    # Forward iteration
    for per=1:Tran_p
        # get distribution
        W = W_seq[per]
        D = reshape(D_seq[per,:],1,NJKT)
        D_new=D*W
        D_seq[per+1,:] = D_new./sum(D_new)   
        diff = maximum(abs.(D_new .- D))
   
        # Housing asset market clearing condition
        h_seq = pol_h_seq[:,per]
        H_d = dot(D_new, h_seq)
        Hmcc_error[per] = H_d - param["H_s"]

        # rental market clearing condition
        s_seq = pol_s_seq[:,per]
        S_d = dot(D_new, s_seq)
        Rmcc_error[per] = S_d - param["H_s"]

        # gini
        tol_cut = 1e-16
        a_seq = pol_a_seq[:,per]
        c_seq = pol_c_seq[:,per]
        q = q_seq[per]
        ρ = ρ_seq[per]
        income = vec(Y_grid) + r*a_seq.*(a_seq.>0) + (ρ-δ_R*q)*(h_seq.-s_seq).*(h_seq.>s_seq)
        w = a_seq + q*h_seq
        gini_seq[per, 1] = Gini(D_seq[per+1,:], income, tol_cut)
        gini_seq[per, 2] = Gini(D_seq[per+1,:], w, tol_cut)
        gini_seq[per, 3] = Gini(D_seq[per+1,:], q*h_seq, tol_cut)
        gini_seq[per, 4] = Gini(D_seq[per+1,:], c_seq, tol_cut)
        

    end

    writedlm("$(Dir)/output/transition/q_seq.txt", q_seq, '\t')
    writedlm("$(Dir)/output/transition/rho_seq.txt", ρ_seq, '\t')
    writedlm("$(Dir)/output/transition/D_seq.txt", D_seq, '\t')
    writedlm("$(Dir)/output/transition/h_seq.txt", pol_h_seq, '\t')
    writedlm("$(Dir)/output/transition/a_seq.txt", pol_a_seq, '\t')
    writedlm("$(Dir)/output/transition/s_seq.txt", pol_s_seq, '\t')
    writedlm("$(Dir)/output/transition/c_seq.txt", pol_c_seq, '\t')
    writedlm("$(Dir)/output/transition/gini_seq.txt", gini_seq, '\t')
    writedlm("$(Dir)/output/transition/V_tr.txt", V_tr, '\t')
    writedlm("$(Dir)/output/transition/Hmcc_error.txt", Hmcc_error, '\t')
    writedlm("$(Dir)/output/transition/Rmcc_error.txt", Rmcc_error, '\t')

    return Hmcc_error, Rmcc_error, gini_seq
end

function moving_average(arr)
    n = length(arr)
    result = zeros(n) # 결과를 저장할 배열
    
    for i in 1:n
        if i == 1
            # 첫 번째 값은 그대로 사용
            result[i] = arr[i]
        elseif i == n
            # 마지막 값은 그대로 사용
            result[i] = arr[i]
        else
            # 가운데 값은 양옆 값을 포함한 평균
            result[i] = (arr[i-1] + arr[i] + arr[i+1]) / 3
        end
    end
    
    return result
end

###############################################################################
# iteration

# given (q, ρ) path
Tran_p = 60
Tr   =  0.137962382

# start with linear guess
q_start, q_end = 1.77904287493328, 1.78135986422181
ρ_start, ρ_end = 0.130186284267784, 0.130264346011243
gini_start = [0.41906809790096256 0.5987413491884841 0.7189663884330421 0.3620223590342342]
#q_ini = log_type_guess(q_start, q_end, Tran_p)
#ρ_ini = log_type_guess(ρ_start, ρ_end, Tran_p)
q_ini = readdlm("$(Dir)/output/transition/old/q_seq.txt", '\t')
ρ_ini = readdlm("$(Dir)/output/transition/old/rho_seq.txt", '\t')
q_ini = moving_average(q_ini)
ρ_ini = moving_average(ρ_ini)
#weight = log_type_guess(1, 0, Tran_p)
weight = range(1, 0.5, Tran_p)
weight = weight.^2

# V_end (SS V with HRET=0)
V_end = readdlm("$(Dir)/output/no_HRET/V.txt", 't')
V_end = reshape(V_end, N, J, K, T)
# initial distribution (SS D with HRET=1)
D_ini = readdlm("$(Dir)/output/benchmark/D.txt", 't')
# D_ini = readdlm("$(Dir)/output/no_HRET/D.txt", 't')

global q_seq = q_ini
global ρ_seq = ρ_ini

println("=============  Start ========================")

for tran_iter=1:20
    Hmcc_error, Rmcc_error, gini_seq = Transition_one(q_seq, ρ_seq, V_end, D_ini, Tran_p, Grid_param)
    max_error = maximum(abs.([Hmcc_error; Rmcc_error]))
    println("=======================================================")
    println("iteration: $(tran_iter), maximum error: $(max_error)")
    println("=======================================================")
    if max_error < 1e-4
        println("converged!")
        break
    end

    # graph
    time = range(0, Tran_p, Tran_p+1)
    q_dyn = [q_start ; q_seq] 
    ρ_dyn = [ρ_start; ρ_seq]
    gini_dyn = [gini_start; gini_seq]
    
    f1 = plot(time, 3*q_dyn, label="", title="House price", titlefontsize=10, color=:blue, linewidth=2)
    plot!(f1, time, 3*q_dyn[1]*ones(Tran_p+1), label="", linestyle=:dash, color=:black)
    plot!(f1, time, 3*q_dyn[end]*ones(Tran_p+1), label="", linestyle=:dash, color=:black)

    f2 = plot(time, ρ_dyn, label="", title="Rent", titlefontsize=10, color=:blue, linewidth=2)
    yaxis!(:log, formatter = x -> string(round(x, digits=4)))
    plot!(f2, time, ρ_dyn[1]*ones(Tran_p+1), label="", linestyle=:dash, color=:black)
    plot!(f2, time, ρ_dyn[end]*ones(Tran_p+1), label="", linestyle=:dash, color=:black)

    f3 = plot(time, gini_dyn[:,1], label="", title="Income Gini", titlefontsize=10, color=:blue, linewidth=2)
    plot!(f3, time, gini_dyn[1,1]*ones(Tran_p+1), label="", linestyle=:dash, color=:black)
    plot!(f3, time, gini_dyn[end,1]*ones(Tran_p+1), label="", linestyle=:dash, color=:black)

    f4 = plot(time, gini_dyn[:,2], label="", title="Wealth Gini", titlefontsize=10, color=:blue, linewidth=2)
    plot!(f4, time, gini_dyn[1,2]*ones(Tran_p+1), label="", linestyle=:dash, color=:black)
    plot!(f4, time, gini_dyn[end,2]*ones(Tran_p+1), label="", linestyle=:dash, color=:black)

    f5 = plot(time, gini_dyn[:,3], label="", title="Housing asset Gini", titlefontsize=10, color=:blue, linewidth=2)
    plot!(f5, time, gini_dyn[1,3]*ones(Tran_p+1), label="", linestyle=:dash, color=:black)
    plot!(f5, time, gini_dyn[end,3]*ones(Tran_p+1), label="", linestyle=:dash, color=:black)

    f6 = plot(time, gini_dyn[:,4], label="", title="Consumption Gini", titlefontsize=10, color=:blue, linewidth=2)
    plot!(f6, time, gini_dyn[1,4]*ones(Tran_p+1), label="", linestyle=:dash, color=:black)
    plot!(f6, time, gini_dyn[end,4]*ones(Tran_p+1), label="", linestyle=:dash, color=:black)

    plot(f1, f2, f3, f4, f5, f6, layout=(2, 3), size=(1000, 600), margin=10Plots.px, legend=:topleft)
    savefig("$(Dir)/output/transition/transition_path_$(tran_iter).png")

    ################### error ################################################
    g = plot(time[2:end], Hmcc_error, label="Housing market error", title="Market clearing error", titlefontsize=10, color=:red, linewidth=2, ylim=(-0.0005, 0.0005))
    plot!(g, time[2:end], Rmcc_error, label="Rental market error", color=:blue, linewidth=2)
    plot!(g, time[2:end], zeros(Tran_p), label="", linestyle=:dash, color=:black)
    plot(g, layout=(1, 1), size=(300, 300), margin=10Plots.px, legend=:topleft)
    savefig("$(Dir)/output/transition/error_$(tran_iter).png")


    # update
    global q_seq = q_seq + 0.20*(weight.*Hmcc_error)
    global ρ_seq = ρ_seq + 0.05*(weight.*Rmcc_error)

    q_ma = moving_average(q_seq)
    ρ_ma = moving_average(ρ_seq)
    global q_seq = q_ma
    global ρ_seq = ρ_ma

end
