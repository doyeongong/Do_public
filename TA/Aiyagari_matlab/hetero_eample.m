% Problem Set 3 Q5
% Yeongwoong Do
clear
clc

tic  % 시간측정 시작

%% Parameters
a_min=0;
a_max=15;
N_a = 100;
beta = 0.95;
r   = 0.035;
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

toc % 시간측정 끝

fprintf("반복횟수 %d회에서 오차 %.6f로 수렴하였습니다.\n", iter, diff)

% mesh grid
[Y_mat, A_mat]=meshgrid(y_grid, a_grid');
pol_c = Y_mat + (1+r)*A_mat - pol_a;

save("value.mat","V_new");
save("policy_a.mat","pol_a");
save("policy_c.mat","pol_c");



%% Graph for (b)

figure(1)
clf

sgtitle("Results of VFI")

% value function
subplot(1,3,1)
plot(a_grid, V_new, 'Linewidth', 2)
title("Value function")
xlabel("a")
ylabel("v(a,y)")
legend("y=low","y=middle","y=high","Location","best")

% policy function (a')
subplot(1,3,2)
plot(a_grid, pol_a, 'Linewidth', 2)
hold on
plot(a_grid, a_grid, '--k', 'Linewidth',1)
hold off
title("Policy function (a')")
xlabel("a")
ylabel("a'(a,y)")
legend("y=low","y=middle","y=high","","Location","best")

% policy function (c)
subplot(1,3,3)
plot(a_grid, pol_c, 'Linewidth', 2)
hold on
plot(a_grid, a_grid, '--k', 'Linewidth',1)
hold off
title("Policy function (c)")
xlabel("a")
ylabel("c(a,y)")
legend("y=low","y=middle","y=high","","Location","best")

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

fprintf(" Average Income (Y) in the Analytical Solution   : %.4f \n", Y_ans)
fprintf(" Average Income (Y) in the Monte Carlo Simulation: %.4f \n", Y_sim)
fprintf(" Average Asset  (A) in the Monte Carlo Simulation: %.4f \n", A_sim)

figure(2)
clf
time = linspace(1,Sim_T,Sim_T);
sgtitle("Aggregation by the Monte Carlo simulation method")

subplot(1,2,1)
plot(time,Y_ser, '-b','Linewidth', 0.5)
hold on
plot(time,Y_ans*ones(Sim_T,1), '-r','Linewidth', 1)
plot(time(burn_in+1:end),Y_sim*ones(Sim_T-burn_in,1), '-c','Linewidth', 1)
xline(burn_in+1,'--k','Linewidth',1)
hold off
title("Y simulation")
legend("Y series","Y analytical mean","Y simulation mean","burn in","Location","best")

subplot(1,2,2)
plot(time,A_ser, '-b','Linewidth', 0.5)
hold on
plot(time(burn_in+1:end),A_sim*ones(Sim_T-burn_in,1), '-c','Linewidth', 1)
xline(burn_in+1,'--k','Linewidth',1)
hold off
title("A simulation")
legend("A series","A simulation mean","burn in","Location","best")













