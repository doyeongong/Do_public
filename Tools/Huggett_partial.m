function [A, result] = Huggett_partial(q)
% Huggett(1993, JEDC) Model
% Yeongwoong Do (2024-04-15)
% step 1~2 in the algorithm of Huggett

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

fprintf('-------------------------------\n');
fprintf('채권 가격(q) : %.4f \n', q);
    
%% (Step1) Get policy function from the VFI

% prepare for the PFI
V_ini = zeros(N,K);             %intial guess for the policy function
V_old = V_ini;
V_new = zeros(N,K);
pol_a_prime = zeros(N,K);

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

% CDF
psi_mat = [psi_new(1:N), psi_new(N+1:end)];
Psi_mat = psi_mat;
for i=2:N
    Psi_mat(i,:)=Psi_mat(i-1,:)+psi_mat(i,:);
end

% aggregation
A = sum(a_grid'*psi_mat);

fprintf('총 채권 수요(A) : %.4f \n', A);

% save results
result = struct();
result.valuef = V_new;
result.policyf = pol_a_prime;
result.stat_den = psi_mat;
result.stat_dist = Psi_mat;

%% sub-function

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

function transition_matrix = BTM(a_grid, a_prime)

% transition matrix 초기화
n = length(a_grid);
transition_matrix = zeros(n, n);

% 각 a_prime 값에 대해 transition matrix 채우기
for i = 1:length(a_prime)
    % a_prime 값이 a_grid 범위 내에 있는지 확인
    if a_prime(i) >= a_grid(1) && a_prime(i) <= a_grid(end)
        % a_prime이 a_grid 내에 있으면 해당하는 인덱스 찾기
        [~, idx] = min(abs(a_grid - a_prime(i)));
        
        % transition matrix 업데이트
        if idx <= n % a_prime 값이 a_grid의 범위 내에 있으면
            % a_prime이 a_grid 값과 정확히 일치할 때
            if a_prime(i) == a_grid(idx)
                transition_matrix(i, idx) = 1;
            elseif a_prime(i) > a_grid(idx) % 가까운 쪽이 왼쪽에 있는 경우
                ratio = (a_prime(i) - a_grid(idx)) / (a_grid(idx+1) - a_grid(idx));
                transition_matrix(i, idx) = 1 - ratio;
                transition_matrix(i, idx+1) = ratio;
            elseif a_prime(i) < a_grid(idx)  % 가까운 쪽이 오른쪽에 있는 경우
                ratio = (a_prime(i) - a_grid(idx-1)) / (a_grid(idx) - a_grid(idx-1));
                transition_matrix(i, idx) = ratio;
                transition_matrix(i, idx-1) = 1-ratio; 
            end
        else % a_prime 값이 a_grid의 최댓값보다 큰 경우
            transition_matrix(i, i) = 1;
        end
    end
end

end

end