% Huggett(1993, JEDC) Model with HACT
% Yeongwoong Do
% 2024-04-18
% Key reference AHLLM(2021) and Code by SeHyoun Ahn
% Note: Moll 코드에서 partial과 GE구하기 파라미터가 바뀜!

clear
clc

fprintf('--------Huggett(1993) by HACT ----------\n');

% 코드 시작 시각 기록
start_time = tic;
start_time_str = datestr(datetime('now'),'yyyy-mm-dd HH:MM:SS');

%%  parameter

sigma = 2;
rho = 0.05; % beta = 1/(1+rho)
e1  = 0.1;
e2  = 0.2;
e_grid   = [e1, e2];
lambda1 = 1.2;
lambda2 = 1.2;
lambda = [lambda1, lambda2];
a_lb = -0.15;
a_ub = 5;
N = 1000;    % number of a_grid (150)
K = 2;      % number of e_grid
a_grid = linspace(a_lb, a_ub, N)';
da = (a_ub - a_lb)/(N-1);
delta = 1000;
% meshgrid 생성 (A_grid = [a_grid, a_grid], E_grid = [e_grid'이 N개 행])
[E_grid, A_grid] = meshgrid(e_grid, a_grid);

% interest rate를 위한 outer loop
r_bracket = [0.01, 0.04];
r = 0.03; % initial guess
max_iter_bisection = 20;
tol_bisection = 1e-5;
diff_bisection = 10;
count_bisection = 0;

while count_bisection < max_iter_bisection && diff_bisection > tol_bisection


%% HJB 풀이

%initial guess(v = v_prime을 의미)
v_ini = zeros(N,K);
% v_ini의 다른 버전: v_ini= u(c)/rho
v_ini = (E_grid+r.*A_grid).^(1-sigma)./(rho*(1-sigma));
v_old = v_ini;
v_new = zeros(N,K);
v_f = zeros(N,K);
v_b = zeros(N,K);

max_iter = 100;
tol = 1e-6;
count = 0;
diff = 10;
A_lambda =[zeros(N,N), diag(lambda1*ones(N,1)); diag(lambda2*ones(N,1)), zeros(N,N)];

while count < max_iter && diff > tol
    
    A = zeros(N*K,N*K);
    % (Step 1) v'.F v'.B 계산

    % forward difference
    v_f(1:N-1,:) = (v_old(2:N,:) - v_old(1:N-1,:))./da;
    v_f(N,:) = (e_grid + r.*a_ub).^(-sigma);
    % backward difference
    v_b(2:N,:) = (v_old(2:N,:) - v_old(1:N-1,:))./da;
    v_b(1,:) = (e_grid + r.*a_lb).^(-sigma);


    % (Step2) c와 s를 계산

    % forward
    c_f = (v_f).^(-1/sigma);
    s_f = E_grid + r.*A_grid - c_f;
    s_f_plus = max(s_f, 0);
    % backward
    c_b = (v_b).^(-1/sigma);
    s_b = E_grid + r.*A_grid - c_b;
    s_b_minus = min(s_b, 0);
    % ss 값
    c_ss = E_grid + r.*A_grid;
    dv_ss = c_ss.^(-sigma);
    Is_f = (s_f>0);
    Is_b = (s_b<0);
    either_not = ones(N,K) - Is_f - Is_b;
    dv = v_f.*(s_f>0) + v_b.*(s_b<0) + dv_ss.*either_not;
    c = dv.^(-1/sigma);

    % (Step3) B v^n+1 = b
    % x, y, z
    X = - s_b_minus./da;
    Y = (- s_f_plus./da + s_b_minus./da)-lambda.*ones(N,1);
    Z = s_f_plus./da;
    % A 행렬 만들기 (spdiags함수를 안쓰는 경우)
    A = diag(Y(:));
    for i=2:N
        A(i,i-1)=X(i,1);
        A(N+i,N+i-1)=X(i,2);
    end
    for i=1:N-1
        A(i,i+1)=Z(i,1);
        A(N+i,N+i+1)=Z(i,2);
    end
    A = A + A_lambda;
    % B 행렬 만들기
    B = (rho+1/delta)*eye(N*K) - A;
    % b 행렬 만들기
    u_mat = c.^(1-sigma)/(1-sigma);
    b_mat = u_mat + (1/delta).*v_old;
    b = b_mat(:);
    % v_new 계산
    v_new_vec = B \ b;

    % (step 4) update
    diff = max(abs(v_old(:) - v_new_vec));
    v_new = reshape(v_new_vec,N,K);
    count = count +1;

    v_old = v_new;

end

%% FPE(KFE)

AT = A';
zero_vec = zeros(N*K,1);
% singular를 피하기 위해 값을 하나주고, 그에 상응 하는 AT 행 변환
zero_vec(1) = 0.1;
AT(1,:) = [1, zeros(1,N*K-1)];
% 안정적 분포 구하기
g = AT \ zero_vec;
% renormalize
g = g./(sum(g)*da);
g_mat = [g(1:N), g(N+1:end)];


%% Upadate

% 총저축
S = g_mat(:,1)'*a_grid*da + g_mat(:,2)'*a_grid*da;
diff_bisection = abs(S-0);

fprintf('총 저축(S) : %.4f \n', S);

if S>0 % excess demand
    r_bracket(2)=r;
elseif S<0 % excess supply
    r_bracket(1)=r;
end
r = sum(r_bracket)/2;
count_bisection = count_bisection + 1;   

end


% 총 경과 시간 계산
elapsed_time = toc(start_time);
% 분과 초로 분리
minutes_elapsed = floor(elapsed_time / 60); % 분 계산
seconds_elapsed = mod(elapsed_time, 60);    % 초 계산
% 결과 출력
fprintf('시작 시각: %s\n', start_time_str);
fprintf('종료 시각: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('총 경과 시간: %d분 %.2f초\n', minutes_elapsed, seconds_elapsed);


%% figure

% 저축함수
s_pol = E_grid + r*A_grid - c;

% 분포함수
G = reshape(g./sum(g),N,K);
cumG = zeros(N,K);
cumG(1,:) = G(1,:);
for i=2:N
    cumG(i,:)=cumG(i-1,:)+G(i,:);
end

a_ub_fig = 1;

% Figure
figure(1);
clf;
plot(a_grid,s_pol(:,1),'-b','Linewidth',2)
hold on
plot(a_grid,s_pol(:,2),'-.r','Linewidth',2)
yline(0)
hold off
title('Saving function')
ylabel('s(a,e)')
xlabel('a')
xlim([a_lb, a_ub_fig])
legend('e=0.1','e=0.2','','','Location','best')


% distribution
figure(2);
clf;
plot(a_grid,cumG(:,1),'-b','Linewidth',2)
hold on
plot(a_grid,cumG(:,2),'-.r','Linewidth',2)
hold off
title('Distribution')
ylabel('mass(share)')
xlabel('a')
xlim([a_lb, a_ub_fig])
legend('e=0.1','e=0.2','Location','best')

% density
figure(3);
clf;
plot(a_grid,G(:,1),'-b','Linewidth',2)
hold on
plot(a_grid,G(:,2),'-.r','Linewidth',2)
hold off
title('Density')
ylabel('density')
xlabel('a')
xlim([a_lb, a_ub_fig])
legend('e=0.1','e=0.2','Location','best')
