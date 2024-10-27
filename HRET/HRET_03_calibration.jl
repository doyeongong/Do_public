# New version of HRET(High-value Real Estate Tax)
# edited by Yeongwoong Do (2024.09.05)
# mu -> a, separable utility
# Objective: calibration

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
        "h_min"   =>  1.20,           # lower bound for h
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
        "adj"     =>  1.100            #################  housing value adjustment for property tax (HRET target) (낮출수록 내는 비율 낮아짐)
);

##############  3. Grid #####################################################
a_grid         = densezero(param["a_min"], param["a_max"], param["N_a"])
h_grid = [0.0, 1.22, 2.0, 3.0, 4.07, 5.825, 8.0, 10.0] 
# h_grid         = [0.0; range(param["h_min"], param["h_max"], param["N_J"]-1)]
# h_grid         = [0.0; log_type_grid(param["h_min"], param["h_max"], param["N_J"]-1)]
println(h_grid)
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

function Housing_PE_LC(q, ρ, Tr, T, Grid_param)   
       
    # parameters
    N, J, K, T, a_min, a_max, unit, Π, s_prob, h_grid, a_grid = Grid_param
    util_param    = param["β"], param["α"], param["σ"], param["ξ"], param["b_bar"] , param["ψ"]
    housing_param = param["γ_s"], param["m"], param["r"], param["θ"], param["κ"], param["δ"], param["r"]+param["κ"] , param["δ_R"], param["ϕ"]
    p, R, adj, minus_infi = param["p"], param["R"], param["adj"], param["-∞"]
    γ_s, m, r, θ, κ, δ, r_m, δ_R, ϕ = housing_param
    
    # empty-set
    V, pol_a, pol_c, pol_s, pol_h, total_tax = zeros(N,J,K,T), zeros(N,J,K,T), zeros(N,J,K,T), zeros(N,J,K,T), zeros(N,J,K,T), zeros(N,J,K,T)
    V_itp = Dict()
    for k=1:K
        for j=1:J
            interp = LinearInterpolation(a_grid, V[:,j,k,end])
            V_itp["func$(j)_$(k)"] = interp
        end
    end

@time begin
            
    @showprogress for t=T:-1:1  # Backward!
            # age-specific housing tax vector
            h_tax = (property_tax.(q*h_grid*unit*adj).+ HRET.(q*h_grid*unit*adj,t))*q.*h_grid

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
                        temp_idx                         = argmax(W_temp)
                        pol_a[a_idx, h_idx, z_idx, t]   = a_temp[temp_idx]
                        pol_c[a_idx, h_idx, z_idx, t]   = c_temp[temp_idx]
                        pol_s[a_idx, h_idx, z_idx, t]   = s_temp[temp_idx]
                        pol_h[a_idx, h_idx, z_idx, t]   = h_temp[temp_idx]                      
                        total_tax[a_idx, h_idx, z_idx, t] = τ_y*tax_base 
                        + τ_R_in*(ρ-δ_R*q)*(h_temp[temp_idx]-s_temp[temp_idx])*(temp_idx==3)
                        + (property_tax(q*h_temp[temp_idx]*unit*adj)+ HRET(q*h_temp[temp_idx]*unit*adj,t))*q*h_temp[temp_idx]
                        - R_de*ρ*s_temp[temp_idx]*(temp_idx==1)
                        - mort_credit(a_temp[temp_idx], q, h_temp[temp_idx], unit, r_m)*(temp_idx==2)
                    end
                end
            end
            
            # 다음 연령에서 사용할 interpolation 만들기
            V_itp = Dict()
            for k=1:K
                for j=1:J
                interp = LinearInterpolation(a_grid, V[:,j,k,t])
                V_itp["func$(j)_$(k)"] = interp
                end
            end
            
        end
        
        # 시간측정 end
    end
    
    return V, pol_a, pol_h, pol_s, pol_c, total_tax
end

#################### 6. Market equilibrium ###############################

# given q(H_s) -> rho

function market_error(ρ, q; print_on=1)
    Tr   =  0.137962382
    V, pol_a, pol_h, pol_s, pol_c, total_tax = Housing_PE_LC(q, ρ, Tr, T, Grid_param);
    D = Aggregation(pol_a, pol_h, Grid_param, stat_dist; print_on=0);

    S_d = dot(vec(pol_s).*vec(pol_h.==0),D)
    S_s = dot(vec(pol_h.-pol_s).*vec(pol_h.>pol_s),D)
    H_d = dot(vec(pol_h),D)
    homeownership = dot(vec(pol_h.>0),D)
    println("x: $(ρ), excess: $(S_d - S_s), H_d: $(H_d)")

    writedlm("$(Dir)/output/V.txt", V, 't')
    writedlm("$(Dir)/output/pol_a.txt", pol_a, 't')
    writedlm("$(Dir)/output/pol_h.txt", pol_h, 't')
    writedlm("$(Dir)/output/pol_s.txt", pol_s, 't')
    writedlm("$(Dir)/output/pol_c.txt", pol_c, 't')
    writedlm("$(Dir)/output/D.txt", D, 't')

            
    if print_on==1
        cutoff = 60000/(param["unit"]*param["adj"]*q)
        momentum(D, pol_h, pol_s, pol_a, pol_c, total_tax, cutoff, H_grid, Y_grid, ρ, q, housing_param, Grid_param, R; tol_cut=1e-16)
        age_profile(D, N, J, K, T, q, pol_a, pol_c, pol_h, pol_s, cutoff, Dir; output_print=0)
        h_dist(D, N, J, K, T, pol_h, h_grid; print_on=1)
    end
    
    return S_d - S_s, H_d, homeownership
    
end

function housing_market_error(q)
    function rental(x)
        excess, _, _ = market_error(x, q; print_on=1)
        return excess
    end
    ρ_sol = newton_one_root_secant(0.13010, 0.13020, rental; tol=1e-4, max_iter=15)
    # ρ_sol = bisection_one_root([0.130154, 0.130162], rental; tol=1e-4, max_iter=15)
    return ρ_sol
end


q_ex = 1/(param["unit"]/10000)
println("q:$(q_ex)")
ρ_sol = housing_market_error(q_ex)
println("q: $(q_ex), ρ: $(ρ_sol)")

