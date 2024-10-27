# calibration module

module fitting

using LinearAlgebra, Plots, Statistics, DelimitedFiles

Dir = "C:/Users/82106/Documents/4. 2024학년도/HRET/Model"
include("$(Dir)/A6_Inequality.jl")
using .Inequality

include("$(Dir)/A5_Korean_tax.jl")
using .Korean_tax

export momentum, age_profile, h_dist

function momentum(D, pol_h, pol_s, pol_a, pol_c, total_tax, cutoff, H_grid, Y_grid, ρ, q, housing_param, Grid_param, R ; tol_cut=1e-10)
    N, J, K, T, a_min, a_max, unit, Π, s_prob, h_grid, a_grid = Grid_param
    γ_s, m, r, θ, κ, δ, r_m, δ_R, ϕ = housing_param

    a_vec, h_vec, c_vec, s_vec = vec(pol_a), vec(pol_h), vec(pol_c), vec(pol_s) 
    
    # age grid
    age = zeros(N, J, K, T)
    for t=1:T
        age[:,:,:,t].=t
    end
    age_vec = vec(age);

    homeownership = dot((h_vec.>0),D)
    landlord      = dot((h_vec.>s_vec),D)
    share_HRET    = dot((h_vec.>cutoff),D)
    trading       = dot((pol_h.!=H_grid),D)
    tax_revenue   = dot(vec(total_tax),D)
    old_mat= zeros(N,J,K,T)
    old_mat[:,:,:,R:end].= 1
    
    income = Y_grid + r*pol_a.*(pol_a.>0) + (ρ-δ_R*q)*(pol_h.-pol_s).*(pol_h.>pol_s)
    e_gini = Gini(D, vec(Y_grid), tol_cut)
    i_gini = Gini(D, vec(income), tol_cut)
    w_gini = Gini(D, (a_vec.+q.*h_vec), tol_cut)
    h_gini = Gini(D, q*h_vec, tol_cut)
    c_gini = Gini(D, c_vec, tol_cut)

    C = dot(c_vec,D)
    S = dot(s_vec,D)
    share_shelter = ρ*S / (C+ρ*S)
    indebted = dot((a_vec.<0),D)
    LTV = (-a_vec)./(q.*h_vec)
    avg_LTV = mean(LTV[a_vec.<0])
    HRET_old = dot(vec(pol_h.>cutoff).*vec(old_mat),D) / share_HRET
    HRET_landlord = dot(vec(pol_h.>cutoff).*(h_vec.>s_vec),D) / share_HRET
    mean_earnings = dot(vec(Y_grid), D)
    mean_rent_pay = dot(s_vec.*(h_vec.==0),D)*ρ/sum(D.*(h_vec.==0))

    prop_re = sum(property_tax.(q*h_vec*unit*1.1))
    HRET_re = sum(HRET.(q*h_vec*unit*1.1, age_vec))
    HRET_prop = HRET_re/prop_re

    println("---------------------------------")
    println("ownership : $(homeownership)")
    println("landlord  : $(landlord)")
    println("HRET      : $(share_HRET)")
    println("trading   : $(trading)")
    println("---------------------------------")
    println("earning gini     : $(e_gini)")
    println("income gini      : $(i_gini)")
    println("wealth gini      : $(w_gini)")
    println("housing gini     : $(h_gini)")
    println("consumption gini : $(c_gini)")
    println("---------------------------------")
    println("share ρ*S     : $(share_shelter)")
    println("HRET 61+      : $(HRET_old)")
    println("HRET landlord : $(HRET_landlord)")
    println("HRET/property : $(HRET_prop)")
    println("---------------------------------")
    println("total tax             : $(tax_revenue)")
    println("mean rent pay         : $(mean_rent_pay)")
    println("mean earings          : $(mean_earnings)")
    println("mean consumption      : $(C)")
    println("fraction of indebted  : $(indebted)")
    println("average LTV           : $(avg_LTV)")
    println("---------------------------------")

    return homeownership, landlord, share_HRET, trading, w_gini, c_gini, share_shelter, HRET_old, HRET_landlord
end


function age_profile(D, N, J, K, T, q, pol_a, pol_c, pol_h, pol_s, cutoff, Dir; output_print=0)
    
    D_mat = reshape(D, N, J, K, T)
    C = dot(vec(D), vec(pol_c))

    # By age
    ownership_age= zeros(T,1)
    landlord_age= zeros(T,1)
    HRET_age= zeros(T,1)
    c_age = zeros(T,1)
    w_age = zeros(T,1)
    agg_c = dot(vec(pol_c),D)
    agg_w = dot(vec(pol_a.+q*pol_h),D)
    
    for t=1:T
        ownership_age[t] = sum((pol_h[:,:,:,t].> 0) .* D_mat[:,:,:,t])/sum(D_mat[:,:,:,t]) 
        landlord_age[t] = sum((pol_h[:,:,:,t].>pol_s[:,:,:,t]).* D_mat[:,:,:,t])/sum(D_mat[:,:,:,t]) 
        HRET_age[t] = sum((pol_h[:,:,:,t].>=cutoff).* D_mat[:,:,:,t])/sum(D_mat[:,:,:,t]) 
        c_age[t] = (sum(pol_c[:,:,:,t].*D_mat[:,:,:,t])./sum(D_mat[:,:,:,t]))/agg_c
        w_age[t] = (sum((pol_a[:,:,:,t].+q*pol_h[:,:,:,t]).*D_mat[:,:,:,t])./sum(D_mat[:,:,:,t]))/agg_w
    end

    # by age grid
    age_grid = [25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80];
    N_age = length(age_grid)-1
    owner_by_age_grid = zeros(N_age)
    landlord_by_age_grid = zeros(N_age)
    HRET_by_age_grid = zeros(N_age)
    c_by_age_grid = zeros(N_age)
    w_by_age_grid = zeros(N_age)

    age = range(25,79,55)

    for i=1:N_age
        idx = (age.>=age_grid[i]).*(age.<age_grid[i+1])
        owner_by_age_grid[i]=mean(ownership_age[idx])
        landlord_by_age_grid[i]=mean(landlord_age[idx])
        HRET_by_age_grid[i]=mean(HRET_age[idx])
        c_by_age_grid[i] = mean(c_age[idx])
        w_by_age_grid[i] = mean(w_age[idx])
    end

    life_profile = [owner_by_age_grid, landlord_by_age_grid, HRET_by_age_grid];
    if output_print==1
        writedlm("Model_life_profile.csv", life_profile, ',')
    end

    # age profile 

    age_housing_data = readdlm("$(Dir)/age_profile/housing.csv", ',');
    age_c_nw_data = readdlm("$(Dir)/age_profile/c_nw.csv", ',');

    age_labels = ["25~29", "30~34", "35~39", "40~44", "45~49", "50~54", "55~59", "60~64", "65~69", "70~74", "75~79"]
    model_style = (linetype=:path, marker=:circle, linestyle=:solid, linewidth=2, markersize=4, markercolor=:red, linecolor=:red)
    data_style = (linetype=:path, marker=:circle, linestyle=:solid, linewidth=2, markersize=4, markercolor=:blue, linecolor=:blue)

    # Consumption volatility 
    #c_adj = std(age_c_nw_data[:,2])/std(c_by_age_grid)

    AH1 = plot(age_labels, owner_by_age_grid; model_style..., label="Model", ylim=(0,1), title="Ownership", titlefontsize=10)
    plot!(AH1, age_labels, age_housing_data[2:end,2]; data_style...,label="Data")
    AH2 = plot(age_labels, landlord_by_age_grid; model_style..., label="", title="Landlord share", titlefontsize=10)
    plot!(AH2, age_labels, age_housing_data[2:end,3]; data_style...,label="")
    AH3 = plot(age_labels, HRET_by_age_grid; model_style..., label="", title="HRET payer share", titlefontsize=10)
    plot!(AH3, age_labels, age_housing_data[2:end,4]; data_style...,label="")
    AH4 = plot(age_labels, c_by_age_grid; model_style..., label="", title="Consumption", titlefontsize=10)
    plot!(AH4, age_labels, age_c_nw_data[:,2]; data_style...,label="")
    AH5 = plot(age_labels, w_by_age_grid; model_style..., label="", title="Wealth", titlefontsize=10)
    plot!(AH5, age_labels, age_c_nw_data[:,3]; data_style...,label="")

    plot(AH1, AH2, AH3, AH4, AH5, layout=(5, 1), size=(650, 1200), margin=10Plots.px, legend=:topleft)
    savefig("$(Dir)/age_profile/model_data.png")

    println("ownership in 75-79: $(owner_by_age_grid[end])")
    
    # percentile
    per_grid = range(0.10, 0.90, 9)
    per_c_nw_data = readdlm("$(Dir)/age_profile/percentile.csv", ',');
    per_c_model = perentile_10(D, vec(pol_c), 0).-C.+*(mean(per_c_nw_data[:,2]))
    per_nw_model = perentile_10(D, vec(pol_a.+q*pol_h), 0)

    per_c = plot(per_grid, per_c_model; model_style..., label="", title="Consumption", titlefontsize=10)
    plot!(per_c, per_grid, per_c_nw_data[:,2]; data_style...,label="")
    per_nw = plot(per_grid, per_nw_model; model_style..., label="", title="Wealth", titlefontsize=10)
    plot!(per_nw, per_grid, per_c_nw_data[:,3]; data_style...,label="")
    plot(per_c, per_nw, layout=(2, 1), size=(700, 400), margin=10Plots.px, legend=:topleft)
    savefig("$(Dir)/age_profile/model_data_percentile.png")

end




function h_dist(D, N, J, K, T, pol_h, h_grid; print_on=1)
    h_dist = zeros(J)
    D_mat = reshape(D, N, J, K, T)
    for i=1:J
        h_dist[i] = sum((pol_h.==h_grid[i]).*D_mat)
    end
    if print_on==1
        println("h : $(h_dist)")
    end
end

end