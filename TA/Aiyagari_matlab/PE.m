function [K_s] = PE(r)
%% Parameters
a_min=0;
a_max=15;
N_a = 100;
beta = 0.95;
z = [-0.5 0 0.5];
P = [0.9 0.08 0.02;
     0.15 0.7 0.15;
     0.05 0.15 0.8];
N_y = length(z);

% while loop 구문을 위한 변수들
iter     = 0;
max_iter = 1000;
diff     = 10;
tol      = 1e-5;

%% VFI(value function iteration)
% a grid
% "linspace(a_min, a_max, N_a)"
a_grid = linspace(a_min,a_max,N_a);
% y grid
y_grid = exp(z);

% Step 1) initial guess
V_ini = zeros(N_a, N_y);
V_old = V_ini;
V_new = V_ini;
pol_a = V_ini;

while (iter < max_iter) && (diff>tol)

% Step2) given V_old, calculate V_new
for a_idx = 1:N_a
    a = a_grid(a_idx);
    
    for y_idx = 1:N_y
        y = y_grid(y_idx);
        RHS_temp = zeros(100,1);
        
        for a_p_idx=1:N_a
            a_p = a_grid(a_p_idx);
            c = y + (1+r)*a - a_p;
            if c < 0
                RHS_temp(a_p_idx) = -1e5;
            else
                EV = 0;
                for k=1:N_y
                    EV = EV + P(y_idx,k)*V_old(a_p_idx,k);
                end
                RHS_temp(a_p_idx) = log(c) + beta*EV;
            end
        end

    % max 계산
    [value, a_p_idx_star] = max(RHS_temp);
    a_p_star = a_grid(a_p_idx_star);
    V_new(a_idx,y_idx) = value;
    pol_a(a_idx,y_idx) = a_p_star;    
    end
end

% Step 3) compare and update
diff = max(abs(V_new-V_old),[],'all');
% update
iter = iter + 1;
V_old = V_new;

end



fprintf("반복횟수 %d회에서 오차 %.6f로 수렴하였습니다.\n", iter, diff)

% mesh grid
[Y_mat, A_mat]=meshgrid(y_grid, a_grid');
pol_c = Y_mat + (1+r)*A_mat - pol_a;

save("value.mat","V_new");
save("policy_a.mat","pol_a");
save("policy_c.mat","pol_c");



%% Aggregation

% Monte Carlo Simulation
Sim_T = 3000;
Sim_N = 1000; 
% 초기값은 (a=0, y=0.6065   a_idx=1, y_idx=1, 로 모두 배정)

psuedo_a_idx_ser = zeros(Sim_N,Sim_T);
psuedo_a_idx_ser(:,1) = 1;

psuedo_y_idx_ser = zeros(Sim_N,Sim_T);
psuedo_y_idx_ser(:,1) = 1;

% second period
[pol_a_idx] = to_idx(pol_a, a_grid);
Cum_P=cumsum(P,2);
rng(20241103);


for t=2:Sim_T
    e_vec = rand(Sim_N,1);
    for i=1:Sim_N
        a_idx = psuedo_a_idx_ser(i,t-1);
        y_idx = psuedo_y_idx_ser(i,t-1);
        psuedo_a_idx_ser(i,t) = pol_a_idx(a_idx,y_idx);
        e = e_vec(i);
        if e<=Cum_P(y_idx,1)
            psuedo_y_idx_ser(i,t) = 1;
        elseif e<=Cum_P(y_idx,2)
            psuedo_y_idx_ser(i,t) = 2;
        elseif e<=Cum_P(y_idx,3)
            psuedo_y_idx_ser(i,t) = 3;
        end
    end
end

%% Aggregation2

[stat_dist] = stationary(P);

% average income
Y_ans = y_grid*stat_dist;
Y_ser = zeros(1,Sim_T);
A_ser = zeros(1,Sim_T);
for t=1:Sim_T
    Y_ser(1,t) = mean(y_grid(psuedo_y_idx_ser(:,t))');
    A_ser(1,t) = mean(a_grid(psuedo_a_idx_ser(:,t))');
end

burn_in = 1000;
Y_sim = mean(Y_ser(burn_in+1:end));
A_sim = mean(A_ser(burn_in+1:end));

fprintf(" Average Asset  (A) in the Monte Carlo Simulation: %.4f \n", A_sim)

K_s = A_sim;

end
