# Module for Rouwenhorst

module Rouwenhorst

export Rouwenhorst_ABRS, Rouwenhorst_KS, stationary  # 외부에서는 Rouwenhorst_ABRS 함수만 사용하도록 설정

# load package
using LinearAlgebra

function trans_mat(rho, N)  
    
    # print error message when N < 2
    if N < 2 
        println("Change N!")
            return nothing
    end
    
    p = (1 + rho) / 2
    # transition matrix
    Theta = [p 1-p; 1-p p]
    if N == 2
        return Theta
    else
        Theta_old = Theta
        for i = 3:N
            # Creating new transition matrix
            Theta_new = p * [Theta_old zeros(i-1, 1); zeros(1, i)] +
                        (1-p) * [zeros(i-1, 1) Theta_old; zeros(1, i)] +
                        (1-p) * [zeros(1, i); Theta_old zeros(i-1, 1)] +
                        p * [zeros(1, i); zeros(i-1, 1) Theta_old]
            row_sums = sum(Theta_new, dims=2) # normalizing by row sums
            Theta_old = Theta_new ./ row_sums
        end
        Theta = Theta_old # Update Theta to the latest Theta_old
    end
    return Theta
end

function stationary(Π)
    D, V = eigen(Π')
    idx = findall(abs.(D.- 1) .< 1e-10)
    if isempty(idx)
        error("There is no eigenvalue=1!")
    end
    dist = V[:, idx[1]]
    stat_dist = dist / sum(dist)
    return stat_dist
end

function Rouwenhorst_ABRS(rho, sigma, N)
    # Based on Auclert, Bardoczy, Rognlie, Staub(2021, Econometrica)

    # Transition matrix = Π
    Π = trans_mat(rho, N)
    
    # stationary distribution for Π = pi
    pi = stationary(Π)

    # set std
    s = range(-1, 1, N)
    mean_s = sum(s)/N
    s = s.*(sigma / sqrt(pi'*((s.-mean_s)).^2))
    
    # normalizing
    e_grid = exp.(s)./(pi'*exp.(s))
    
    return e_grid, Π, pi    
end

function Rouwenhorst_KS(rho, sigma, N)
    # Based on Kopecky and Suen(2010, RED)

    # print error message when N < 2
    if N < 2
        println("Change N!")
        return nothing
    end
    
    p = (1 + rho) / 2
    sigma_z = sigma / sqrt(1 - rho^2)
    psi = sqrt(N - 1) * sigma_z
    
    # state space
    y_1 = -psi
    y_N = psi
    y = range(y_1, stop=y_N, length=N)
    
    # transition matrix
    Theta = [p 1-p; 1-p p]
    if N == 2
        return Theta, y
    else
        Theta_old = Theta
        for i = 3:N
            # Creating new transition matrix
            Theta_new = p * [Theta_old zeros(i-1, 1); zeros(1, i)] +
                        (1-p) * [zeros(i-1, 1) Theta_old; zeros(1, i)] +
                        (1-p) * [zeros(1, i); Theta_old zeros(i-1, 1)] +
                        p * [zeros(1, i); zeros(i-1, 1) Theta_old]
            row_sums = sum(Theta_new, dims=2) # normalizing by row sums
            Theta_old = Theta_new ./ row_sums
        end
        Theta = Theta_old # Update Theta to the latest Theta_old
    end

    pi = stationary(Theta)
    return y, Theta, pi
end


end