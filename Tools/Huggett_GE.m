% Replication of Huggett 

% This m.file needs "Huggett_partial.m" !!!!
% (step3)

clear
clc

%% parameter
e_grid = [1.0; 0.1];
beta = 0.99322;
sigma = 1.5;
Pi = [0.925,    1-0.925;
      1-0.5,    0.5];
a_lb = -2;
a_ub = 4;
N = 150;    % number of a_grid
K = 2;      % number of e_grid
a_grid = linspace(a_lb, a_ub, N)';

% price를 위한 outer loop
q_braket = [0.98, 1.02];
% 옵션 설정
options = optimset('TolX', 1e-5); % 허용 오차를 1e-5으로 설정 (원하는 값으로 조정)
% Huggett_partial 함수 정의
fun = @(q) Huggett_partial(q);
% fzero 함수를 사용하여 A가 0이 되는 q 찾기 (초기 추정값으로 0 사용)
q_sol = fzero(fun, q_braket, options);
    
% reproduce results
[A_sol, result] = Huggett_partial(q_sol);
V_sol = result.valuef;
pol_a_prime = result.policyf;
psi_mat = result.stat_den;
Psi_mat = result.stat_dist;


%% Figure

% Figure
figure(1);
clf;
plot(a_grid,a_grid,'--k')
hold on
plot(a_grid,pol_a_prime(:,1),'-b','Linewidth',2)
plot(a_grid,pol_a_prime(:,2),'-r','Linewidth',2)
hold off
title('Policy function')
ylabel('a_{next}(a,e)')
xlabel('a')
legend('','e=1.0','e=0.1','Location','best')

% distribution
figure(2);
clf;
plot(a_grid,Psi_mat(:,1),'-b','Linewidth',2)
hold on
plot(a_grid,Psi_mat(:,2),'-r','Linewidth',2)
hold off
title('Distribution')
ylabel('mass(share)')
xlabel('a')
legend('e=1.0','e=0.1','Location','best')

% density
figure(3);
clf;
plot(a_grid,psi_mat(:,1),'-b','Linewidth',2)
hold on
plot(a_grid,psi_mat(:,2),'-r','Linewidth',2)
hold off
title('Density')
ylabel('density')
xlabel('a')
legend('e=1.0','e=0.1','Location','best')