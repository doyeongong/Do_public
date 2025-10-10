% Problem Set 3 Q5
% Yeongwoong Do
clear
clc

tic % 시간측정 시작
Num_grid= 20;
r_dom = linspace(0.01,0.07,Num_grid);
K_s_vec = zeros(1,Num_grid);
for n=1:Num_grid
    K_s_vec(n)=PE(r_dom(n));
end
toc % 시간측정 끝

%% capital demand
delta = 0.05;
alpha = 0.33;
TFP   = 1;
K_d_fun = @(r) ((r+delta)/(TFP*alpha))^(1/(alpha-1));
K_d_vec = zeros(1,Num_grid);
for n=1:Num_grid
    K_d_vec(n)=K_d_fun(r_dom(n));
end


%% Graph
figure(1)
clf

plot(K_s_vec, r_dom,'--ob','Linewidth',1)
hold on
plot(K_d_vec, r_dom,'--or','Linewidth',1)
xlabel("K")
ylabel("r")
title("capital market")
legend("K^s (capital supply)","K^d (capitl demand)","Location","best")

%% root finding (MCC)


MCC_error = @(r) K_d_fun(r)-PE(r);

r_sol = fzero(MCC_error,0.04);

fprintf("균형이자율: %.4f \n", r_sol) 
fprintf("균형자본량: %.4f \n", K_d_fun(r_sol)) 


%%
figure(1)
clf
beta = 0.95;
K_s_eq = K_d_fun(r_sol);
plot(K_s_vec, r_dom,'--ob','Linewidth',1)
hold on
plot(K_d_vec, r_dom,'--or','Linewidth',1)
scatter(K_s_eq, r_sol, 100, 'r','filled');
text(K_s_eq+0.5, r_sol, sprintf('equilibrium: (%.4f, %.4f)', K_s_eq, r_sol));
yline(1/0.95-1,'--k',"Linewidth",2)
text(0.5, 0.055,sprintf('RA: (%.4f, %.4f)', K_d_fun(1/beta-1), 1/beta-1));
xlabel("K")
ylabel("r")
title("capital market")
legend("K^s (capital supply)","K^d (capitl demand)","RA(1/beta-1)","Location","southeast")










