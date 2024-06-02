% RBC model IRF using ABRS(2021, Econometrica)
% Yeongwoong Do(2024)
clear
clc

%% Set up

% parameter
% eis, beta, vphi, frisch, delta, alpha
eis = 1;
frisch = 1;
delta = 0.025;
alpha = 0.11;
L_ss = 1;
Y_ss = 1;
r_ss = 0.01;

% steady state
% Note : some parameters are endogeneously calibrated (beta, vphi)
beta = 1/(1+r_ss);
w_ss = (1-alpha);
K_ss = alpha*Y_ss/(r_ss+delta);
I_ss = delta*K_ss;
C_ss = Y_ss - I_ss;
Z_ss = Y_ss/(K_ss^alpha*L_ss^(1-alpha));
vphi = w_ss/(L_ss^frisch*C_ss^eis);

% sequence of Z
% dZ_seq(:,1) <= impact period=1
% dZ_seq(:,2) <= impact period=10 but we know at t=1
T = 50;
News = 10;
time = linspace(1,T,T);
rho_z = 0.8;
z_shock = 0.01;
dZ_seq = zeros(T,2);
Z_seq = zeros(T,2);
dZ_seq(1,1) = z_shock*Z_ss;
Z_seq(1,:) = Z_ss + dZ_seq(1,:);
for t=2:T
    if t ~= News
        dZ_seq(t,:) = rho_z*dZ_seq(t-1,:);
    elseif t==News
        dZ_seq(t,1) = rho_z*dZ_seq(t-1);
        dZ_seq(News,2) = z_shock*Z_ss;
    end
    Z_seq(t,:) = Z_ss + dZ_seq(t,:);
end

%% Jacobian

% simple Jacobian
H_U = zeros(T*2,T*2);
H_Z = zeros(T*2,T);
H_ss = H(K_ss, L_ss, Z_ss, L_ss, Z_ss, K_ss, vphi);
h= 1e-5;

% partial
partial_K_t = (H(K_ss+h, L_ss, Z_ss, L_ss, Z_ss, K_ss, vphi) - H_ss)/h ;
partial_L_t = (H(K_ss, L_ss+h, Z_ss, L_ss, Z_ss, K_ss, vphi) - H_ss)/h ;
partial_Z_t = (H(K_ss, L_ss, Z_ss+h, L_ss, Z_ss, K_ss, vphi) - H_ss)/h ;
partial_L_t1 = (H(K_ss, L_ss, Z_ss, L_ss+h, Z_ss, K_ss, vphi) - H_ss)/h ;
partial_Z_t1 = (H(K_ss, L_ss, Z_ss, L_ss, Z_ss+h, K_ss, vphi) - H_ss)/h ;
partial_K_t_1 = (H(K_ss, L_ss, Z_ss, L_ss, Z_ss, K_ss+h, vphi) - H_ss)/h ;


% Full description of H matrix

for i=1:T
    
    % current = diagonal
    % K_t -> EE, MCC
    H_U(i,i) = partial_K_t(1);
    H_U(T+i,i) = partial_K_t(2);  
    % L_t -> EE, MCC
    H_U(i,T+i) = partial_L_t(1);
    H_U(T+i,T+i) = partial_L_t(2); 
    % Z_t -> EE, MCC
    H_Z(i,i) = partial_Z_t(1);
    H_Z(T+i,i) = partial_Z_t(2);
    
    % lead
    if i>1
        % L_t1 -> EE, MCC
        H_U(i-1,T+i) = partial_L_t1(1);
        H_U(T+i-1,T+i) = partial_L_t1(2);
        % Z_t1 -> EE, MCC
        H_Z(i-1,i) = partial_Z_t1(1);
        H_Z(T+i-1,i) = partial_Z_t1(2);
    end
    
    % lag
    if i<T
        % K_t-1 -> EE, MCC
        H_U(i+1,i) = partial_K_t_1(1);
        H_U(T+i+1,i) = partial_K_t_1(2);
    end
    
end

G = -inv(H_U)*H_Z;

% IRF
dK_seq = G(1:T,1:T)*dZ_seq(:,1)./K_ss*100;
dL_seq = G(T+1:2*T,1:T)*dZ_seq(:,1)./L_ss*100;

dK_seq_news = G(1:T,1:T)*dZ_seq(:,2)./K_ss*100;
dL_seq_news = G(T+1:2*T,1:T)*dZ_seq(:,2)./L_ss*100;


%% Figure
figure(1)
clf;

% conventional shock

subplot(2,3,1)
plot(time, 100*dZ_seq(:,1)./Z_ss, "-b", Linewidth=2)
ylabel("% deviation from ss")
xlabel("quaters")
title("TFP shock")

subplot(2,3,2)
plot(time, dK_seq, "-b", Linewidth=2)
hold on
plot(time, zeros(T,1), "--k")
hold off
ylabel("% deviation from ss")
xlabel("quaters")
title("Capital")

subplot(2,3,3)
plot(time, dL_seq, "-b", Linewidth=2)
hold on
plot(time, zeros(T,1), "--k")
hold off
ylabel("% deviation from ss")
xlabel("quaters")
title("Labor")

% News shock

subplot(2,3,4)
plot(time, 100*dZ_seq(:,2)./Z_ss, "-r", Linewidth=2)
ylabel("% deviation from ss")
xlabel("quaters")
title("TFP shock")

subplot(2,3,5)
plot(time, dK_seq_news, "-r", Linewidth=2)
hold on
plot(time, zeros(T,1), "--k")
hold off
ylabel("% deviation from ss")
xlabel("quaters")
title("Capital")

subplot(2,3,6)
plot(time, dL_seq_news, "-r", Linewidth=2)
hold on
plot(time, zeros(T,1), "--k")
hold off
ylabel("% deviation from ss")
xlabel("quaters")
title("Labor")

%% Subfunction

function H_vec = H(Kt, Lt, Zt, Lt1, Zt1, Kt_, vphi)
    % parameters
    eis = 1; frisch = 1; delta = 0.025; alpha = 0.11; beta=1/(1+0.01);     
    Yt = Zt*Kt_^alpha*Lt^(1-alpha);
    It = Kt - (1-delta)*Kt_;
    rt = alpha*Zt*(Kt_/Lt)^(alpha-1)-delta;
    rt1 = alpha*Zt1*(Kt/Lt1)^(alpha-1)-delta;
    wt = (1-alpha)*Zt*(Kt_/Lt)^alpha;
    wt1 = (1-alpha)*Zt1*(Kt/Lt1)^alpha;
    Ct = (wt/(vphi*Lt^frisch))^(1/eis);
    Ct1 = (wt1/(vphi*Lt1^frisch))^(1/eis);
    H_vec = [Ct^(-eis) - beta*(1+rt1)*(Ct1)^(-eis); Yt - Ct - It];
end







