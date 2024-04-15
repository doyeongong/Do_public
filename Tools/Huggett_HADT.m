% Huggett(1993, JEDC) Model
% Yeongwoong Do
% 2024-04-15
% with bisection method

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

% price를 위한 outer loop
q_braket = [0.98, 1.02];
q = sum(q_braket)/2;
max_iter_bisection = 50;
tol_bisection = 1e-5;
diff_bisection = 10;
count_bisection = 0;

while count_bisection < max_iter_bisection && diff_bisection > tol_bisection

fprintf('-------------------------------\n');
fprintf('채권 가격(q) : %.4f \n', q);
    
%% (Step1) Get policy function from the VFI

% prepare for the PFI
a_grid = linspace(a_lb, a_ub, N)';
V_ini = zeros(N,K);             %intial guess for the policy function
V_old = zeros(N,K);
V_new = zeros(N,K);
pol_a_prime = zeros(N,K);
V_old = V_ini;

% loop
max_iter = 5000;
count = 0;
tol = 1e-4;
diff = 10;

tic;

while count < max_iter && diff > tol
for e_idx=1:K
    for a_idx=1:N
        a = a_grid(a_idx);
        e = e_grid(e_idx);
        f_obj = @(x) minus_RHS(x, a, e, q, a_grid, V_old, Pi, e_idx, beta);
        % a_lb <= a' <= (a+e)/q (C가 음수가 안될 조건)
        a_ub_endo = (a+e)/q;
        [pol_a_prime(a_idx,e_idx), minus_V_new]= fminbnd(f_obj, a_lb, a_ub_endo);
        % V와 a' 수정
        V_new(a_idx,e_idx) = -minus_V_new;
        if minus_V_new == 1e4
            pol_a_prime(a_idx,e_idx)=a_lb;
        end
    end
end
diff = max(abs(V_new(:)-V_old(:)));
count = count + 1;
V_old = V_new;

end

% 루프가 종료되었을 때의 경과 시간 출력
elapsed_time = toc;
fprintf('VFI 반복횟수: %d 회\n', count);
fprintf('VFI 경과시간: %.2f 초\n', elapsed_time);

%% (Step2) Get stationary distribution

% psi_0 (initial guess for the stationary distribution)
psi_ini = zeros(N*K,1);
psi_ini(1) = 1;

a_h_tr = BTM(a_grid,pol_a_prime(:,1));
a_l_tr = BTM(a_grid,pol_a_prime(:,2));
W = [Pi(1,1)*a_h_tr, Pi(1,2)*a_h_tr;
    Pi(2,1)*a_l_tr, Pi(2,2)*a_l_tr];

step2_count = 0;
step2_diff = 10;
step2_tol = 1e-6;
psi_old = psi_ini;
psi_new = zeros(N*K,1);

while step2_count < max_iter && step2_diff>step2_tol
    psi_new = W'*psi_old;  % Notation 주의!  
    step2_diff = max(abs(psi_new-psi_old));
    step2_count= step2_count + 1;
    psi_old = psi_new;
end

% check
% e_stat = eig_stat_dist(Pi);
% e_stat_BTM = [sum(psi_new(1:N)); sum(psi_new(N+1:end))];
% abs(e_stat - e_stat_BTM)

% CDF
psi_mat = [psi_new(1:N), psi_new(N+1:end)];
Psi_mat = psi_mat;
for i=2:N
    Psi_mat(i,:)=Psi_mat(i-1,:)+psi_mat(i,:);
end


%% (Step 3) Update

% aggregation
A = sum(a_grid'*psi_mat);
diff_bisection = abs(A-0);

fprintf('총 채권 수요(A) : %.4f \n', A);

if A>0 % excess demand
    q_braket(1)=q;
elseif A<0 % excess supply
    q_braket(2)=q;
end
q = sum(q_braket)/2;
count_bisection = count_bisection + 1;   


end

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



%% function

function utility = U(c)
    sigma = 1.5;
    utility = c^(1-sigma)/(1-sigma);
end

function error = minus_RHS(x, a, e, q, a_grid, V_old, Pi, e_idx, beta)
    sigma= 1.5;
    U = @(c_arg) c_arg^(1-sigma)/(1-sigma);
    c = a+e-x*q;
    if c <=0
        error = 1e4;
    else
        error = -U(a+e-x*q) - beta*(interp1(a_grid, V_old(:,1), x, 'linear')*Pi(e_idx,1) + interp1(a_grid, V_old(:,2), x, 'linear')*Pi(e_idx,2));
    end
end