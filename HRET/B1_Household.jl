# Module for HA with housimg model help functions
# Yeongwoong Do(2024)

module Household

export U, bequest, EV, Γ, intra_FOC, mort_credit, Renter_problem, Owner_problem, Landlord_problem, Aggregation, tran_Lambda

# -infinity
minus_infi = -1e10

# package
using Optim, SparseArrays
Dir = "C:/Users/82106/Documents/4. 2024학년도/HRET/Model"

include("$(Dir)/A1_Rouwenhorst.jl")
include("$(Dir)/A2_Grid.jl")
include("$(Dir)/A3_transition_matrix.jl")
include("$(Dir)/A4_root_finding.jl")
include("$(Dir)/A5_Korean_tax.jl")
include("$(Dir)/A6_Inequality.jl")

using .Rouwenhorst
using .Grid
using .transition_matrix
using .root_finding
using .Korean_tax
using .Inequality


function U(c, s, η, α, ξ)
    # c: non-housing consumption,  s: shelter service
    return c^(1-η)/(1-η) + α*s^(1-ξ)/(1-ξ)
end

function WY(c, s, σ, η, α, ξ)
    # c: non-housing consumption,  s: shelter service
    u = (c^(1-η) + α*((1-η)/(1-ξ))*s^(1-ξ))^(1/(1-η))
    v = u^(1-σ)/(1-σ)
    return v
end

function bequest(a, B, ψ)
    # a : Net wealth
    return B*(a)^(1-ψ)/(1-ψ)
end

function EV(a_p, j_p, e_idx, V_itp, Π, K)
    # V_itp : Dict
    exp=0
    for k=1:K
        itp = V_itp["func$(j_p)_$(k)"]
        exp += Π[e_idx, k]*itp(a_p)
    end
    return exp
end

# transaction cost
function Γ(h, h_p, q, τ_b, γ_s, δ, m)
    # τ_b : acquisition tax rate
    # γ_s : selling cost
    # m   : moving cost
    if h!=h_p
        cost = (τ_b+m)*q*h_p + γ_s*q*h*(1-δ) 
    else
        cost = 0
    end
   return cost
end

# mortgage deduction
function mort_credit(a_p, q, h_p, unit, r_m)
    debt = (a_p)*(a_p<0)*unit
    house_value = q*h_p*unit
    mort_tr = (mortgage_lumpsum(house_value, debt, r_m))/unit
    return mort_tr
end

function intra_FOC(x, ρ, α, σ, ξ)
    function error(s)
        if s<=0
            return minus_infi
        else
            return (ρ/α)^(1/σ)*s^(ξ/σ)+ρ*s - x
        end
    end
    s_sol = newton_one_root_secant(1e-6, 2e-6, error; tol=1e-8, max_iter=20)
    c_sol = (ρ/α)^(1/σ)*s_sol^(ξ/σ)
    return c_sol, s_sol
end

function Renter_problem(q, ρ, h, t, z_idx, common_x, V_itp, R_de, util_param, Grid_param, housing_param)
    
    γ_s, m, r, θ, κ, δ, r_m, δ_R, ϕ = housing_param
    _, _, K, T, _, a_max, _, Π, s_prob, _, _ = Grid_param
    β, α, σ, ξ, b_bar, ψ = util_param
    bar_x = common_x - Γ(h, 0, q, 0, γ_s, δ, 0) 
    a_upper = minimum([bar_x, a_max])
    if a_upper <= 0.0 # No solution case
        W, pol_a, pol_c, pol_s = minus_infi, 0, 0, 0
    else
        function minus_RHS_rent(a_p)
            c, s = intra_FOC(bar_x-a_p, ρ*(1-R_de), α, σ, ξ)
            if c<=0 || s<=0
                return -(minus_infi)
            else
                u = U(c, s, σ, α, ξ)
                if t==T  
                    return -(u + bequest(a_p, b_bar, ψ))
                else
                    return -(u + β*s_prob[t]*EV(a_p, 1, z_idx, V_itp, Π, K))
                end
            end
        end
        results = Optim.optimize(minus_RHS_rent, 0.0, a_upper)
        W = -Optim.minimum(results)
        if W == minus_infi
            pol_a, pol_c, pol_s = 0, 0, 0
        else
            pol_a = Optim.minimizer(results)
            pol_c, pol_s = intra_FOC(bar_x-pol_a, ρ*(1-R_de), α, σ, ξ)
        end
    end
    return W, pol_a, pol_c, pol_s
end

# owner given j_p
function owner_under_j_p(j_p, h_p, q, h, t, z_idx, bar_x, V_itp, DTI_lim, util_param, Grid_param, housing_param)
    γ_s, m, r, θ, κ, δ, r_m, δ_R, ϕ = housing_param
    _, _, K, T, a_min, a_max, unit, Π, s_prob, _, _ = Grid_param
    β, α, σ, ξ, b_bar, ψ = util_param
    a_lower = maximum([-θ*q*h_p, DTI_lim, a_min])
    a_upper = minimum([bar_x-q*h_p, a_max])
    
    if a_upper <= a_lower # No solution case
        W, pol_a, pol_c = minus_infi, 0, 0
    else
        function minus_RHS_owner(a_p)
            c = bar_x - a_p - q*h_p  + mort_credit(a_p, q, h_p, unit, r_m)
            if c<=0
                return -(minus_infi)
            else
                u = U(c, h_p, σ, α, ξ)
                if t==T
                    return -(u + bequest(a_p+q*h_p, b_bar, ψ))
                else
                    return -(u + β*s_prob[t]*EV(a_p, j_p, z_idx, V_itp, Π, K))
                end
            end
        end
        results = Optim.optimize(minus_RHS_owner, a_lower, a_upper)
        W = -Optim.minimum(results)
        if W == minus_infi
            pol_a, pol_c = 0, 0
        else
            pol_a = Optim.minimizer(results)
            pol_c = bar_x - pol_a - q*h_p + mort_credit(pol_a, q, h_p, unit, r_m)
        end
    end
    return W, pol_a, pol_c
end

# j_p loop  (fast version)
function Owner_problem(j_p_lower, q, h, t, z_idx, common_x, V_itp, h_tax, DTI_lim, util_param, Grid_param, housing_param)
    γ_s, m, r, θ, κ, δ, r_m, δ_R, ϕ = housing_param
    _, J, K, T, a_min, a_max, unit, Π, s_prob, h_grid, _ = Grid_param
    v_temp    = zeros(J-(j_p_lower-1))
    a_p_temp  = zeros(J-(j_p_lower-1))
    c_temp    = zeros(J-(j_p_lower-1))
        
    for j_p=j_p_lower:J
        h_p = h_grid[j_p]
        τ_b = acquisition_tax(q*h_p*unit)
        bar_x = common_x - Γ(h, h_p, q, τ_b, γ_s, δ, m) - h_tax[j_p]
        v_temp[j_p-(j_p_lower-1)], a_p_temp[j_p-(j_p_lower-1)], c_temp[j_p-(j_p_lower-1)] = 
        owner_under_j_p(j_p, h_p, q, h, t, z_idx, bar_x, V_itp, DTI_lim, util_param, Grid_param, housing_param)
    end
    
    V = maximum(v_temp)
    
    if V == minus_infi
        pol_h_p, pol_a_p, pol_c = 0, 0, 0
    else
        pol_h_p   = h_grid[argmax(v_temp) + (j_p_lower-1)]
        pol_a_p       = a_p_temp[argmax(v_temp)] 
        pol_c         = c_temp[argmax(v_temp)] 
    end
    
    return V, pol_a_p, pol_c, pol_h_p
end

# landlord given j_p
function landlord_under_j_p(j_p, h_p, q, ρ, h, t, z_idx, bar_x, V_itp, τ_R, DTI_lim, util_param, Grid_param, housing_param)
    γ_s, m, r, θ, κ, δ, r_m, δ_R, ϕ = housing_param
    _, _, K, T, a_min, a_max, unit, Π, s_prob, _, _ = Grid_param
    β, α, σ, ξ, b_bar, ψ = util_param
    a_lower = maximum([-θ*q*h_p, DTI_lim, a_min])
    a_upper = minimum([bar_x-q*h_p, a_max])
    landlord_ρ = (1-τ_R)*(ρ-δ_R*q)

    if a_upper <= a_lower # No solution case
        W, pol_a, pol_c, pol_s = minus_infi, 0, 0, 0
    else
        function minus_RHS_landlord(a_p)
            landlord_x = bar_x - a_p - q*h_p + landlord_ρ*h_p
            c, s = intra_FOC(landlord_x, landlord_ρ, α, σ, ξ)
            if c<=0 || s<=0
                return -(minus_infi)
            else
                u = U(c, s, σ, α, ξ)
                if t==T   
                    return -(u + bequest(a_p+q*h_p, b_bar, ψ))
                else
                    return -(u + β*s_prob[t]*EV(a_p, j_p, z_idx, V_itp, Π, K))
                end
            end
        end
        results = Optim.optimize(minus_RHS_landlord, a_lower, a_upper)
        W = -Optim.minimum(results)
        if W == minus_infi
            pol_a, pol_c, pol_s = 0, 0, 0
        else
            pol_a = Optim.minimizer(results)
            pol_c, pol_s = intra_FOC(bar_x - pol_a - q*h_p + landlord_ρ*h_p, landlord_ρ, α, σ, ξ)
        end
    end
    return W, pol_a, pol_c, pol_s
end

# j_p loop  (fast version)
function Landlord_problem(j_p_lower, q, ρ, h, t, z_idx, common_x, V_itp, h_tax, τ_R, DTI_lim, util_param, Grid_param, housing_param)
    γ_s, m, r, θ, κ, δ, r_m, δ_R, ϕ = housing_param
    _, J, K, T, a_min, a_max, unit, Π, s_prob, h_grid, _ = Grid_param
    v_temp    = zeros(J-(j_p_lower-1))
    a_p_temp  = zeros(J-(j_p_lower-1))
    c_temp    = zeros(J-(j_p_lower-1))
    s_temp    = zeros(J-(j_p_lower-1))
        
    for j_p=j_p_lower:J
        h_p = h_grid[j_p]
        τ_b = acquisition_tax(q*h_p*unit)
        bar_x = common_x - Γ(h, h_p, q, τ_b, γ_s, δ, m) - h_tax[j_p] - ϕ
        v_temp[j_p-(j_p_lower-1)], a_p_temp[j_p-(j_p_lower-1)], c_temp[j_p-(j_p_lower-1)], s_temp[j_p-(j_p_lower-1)] = 
        landlord_under_j_p(j_p, h_p, q, ρ, h, t, z_idx, bar_x, V_itp, τ_R, DTI_lim, util_param, Grid_param, housing_param)
    end
    
    V = maximum(v_temp)
    
    if V == minus_infi
        pol_h_p, pol_a_p, pol_c, pol_s = 0, 0, 0, 0
    else
        pol_h_p   = h_grid[argmax(v_temp) + (j_p_lower-1)]
        pol_a_p   = a_p_temp[argmax(v_temp)] 
        pol_c     = c_temp[argmax(v_temp)] 
        pol_s     = s_temp[argmax(v_temp)] 
    end
    
    return V, pol_a_p, pol_c, pol_s, pol_h_p
end

function Aggregation(pol_a, pol_h, Grid_param, stat_dist; print_on=1)
    # population growth rate
    n = 0
    
    # parameters
    N, J, K, T, _, _, _, Π, s_prob, h_grid, a_grid = Grid_param
    NJK = N * J * K
    NJKT = N * J * K * T

    # Initialize sparse matrix W
    W = spzeros(NJKT, NJKT)
    
    # Statistical matrix
    stat_mat = repeat(stat_dist', outer=(K, 1))

    # 각 블록을 적절한 위치에 배치
    for t in 1:T
        rows = ((t-1)*NJK+1):(t*NJK)
        
        # 죽음: 물려받음 + 생산성 랜덤
        # W[rows, 1:NJK] = BTM_housing(pol_a[:,:,:,t], pol_h[:,:,:,t], stat_mat, N, J, K, a_grid, h_grid) * (1 - s_prob[t] / (1 + n))  
        W[rows, 1:NJK] = BTM_housing(zeros(N,J,K), zeros(N,J,K), stat_mat, N, J, K, a_grid, h_grid).*(1-s_prob[t]/(1+n))      # a=0, h=0에서 시작

        # 생존: 다음 연령으로 전이
        if t < T
            W[rows, (t*NJK+1):(t*NJK+NJK)] = BTM_housing(pol_a[:,:,:,t], pol_h[:,:,:,t], Π, N, J, K, a_grid, h_grid) * (s_prob[t] / (1 + n))
        end
    end

    # psi 초기값
    psi_ini = ones(1, NJKT) ./ NJKT
    
    max_iter = 10000
    tol = 1e-16
    psi_old = psi_ini
    diff = 10
        
    for i in 1:max_iter
        psi_new = psi_old * W
        psi_new ./= sum(psi_new)
        diff = maximum(abs.(psi_new .- psi_old))
        psi_old = copy(psi_new)
        if diff < tol
            if print_on ==1
                println("iteration: $i,  error: $diff")
            end
            break
        end
        if i==max_iter && diff>tol
            println("fail to converge in finidng the stationary distribution")
        end
    end
    
    return vec(psi_old)
end

function tran_Lambda(pol_a, pol_h, Grid_param, stat_dist)
    # population growth rate
    n = 0
    
    # parameters
    N, J, K, T, _, _, _, Π, s_prob, h_grid, a_grid = Grid_param
    NJK = N * J * K
    NJKT = N * J * K * T

    # Initialize sparse matrix W
    W = spzeros(NJKT, NJKT)
    
    # Statistical matrix
    stat_mat = repeat(stat_dist', outer=(K, 1))

    # 각 블록을 적절한 위치에 배치
    for t in 1:T
        rows = ((t-1)*NJK+1):(t*NJK)
        
        # 죽음: 물려받음 + 생산성 랜덤
        # W[rows, 1:NJK] = BTM_housing(pol_a[:,:,:,t], pol_h[:,:,:,t], stat_mat, N, J, K, a_grid, h_grid) * (1 - s_prob[t] / (1 + n))  
        W[rows, 1:NJK] = BTM_housing(zeros(N,J,K), zeros(N,J,K), stat_mat, N, J, K, a_grid, h_grid).*(1-s_prob[t]/(1+n))      # a=0, h=0에서 시작

        # 생존: 다음 연령으로 전이
        if t < T
            W[rows, (t*NJK+1):(t*NJK+NJK)] = BTM_housing(pol_a[:,:,:,t], pol_h[:,:,:,t], Π, N, J, K, a_grid, h_grid) * (s_prob[t] / (1 + n))
        end
    end
    
    return W
end


end