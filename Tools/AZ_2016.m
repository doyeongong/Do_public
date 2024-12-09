%Alpanda and Zubairy(2016, JMCB)
clear
clc

%% Parameters

% 내생변수인것 같지만 실제로 데이터를 사용한것
r = 0.01;          % 분기 1%   (연율 4% 가정)
beta_I = 1/1.015;  % 분기 1.5% (연율 200bp 스프레드 가정)

% exogenous calibration
l_P   = 1;            % normalization
alpha = 0.24;         % capital income share
theta = 1;            % inverse Frisch
rho_b = 0.85;         % persistency of debt
phi   = 0.70;         % average LTV
chi_P = 0.038; chi_I = 0.035; chi_R = 0.015; % transfer share
rrho_b = 0.005;       % transfer response coefficient

% tax rate
tau_b = 0.15; tau_y = 0.3; tau_R = 0.2; 
tau_k = 0.4; tau_p = 0.014/4; tau_c = 0.05; 

%(NIPA(Q), National BS(Y) ratio
ih_y  = 0.05; ik_y = 0.12; g_y  = 0.18;  
h_y   = 5.09;  k_y = 6.00; b_h  = 0.30; hR_h  = 0.2;   


%% Stationary Equilibrium

% (6) <Patient> Government bonds/lending:
beta_P = 1/(1+(1-tau_b)*r);

% (10) <Impatient> Borrowing:
mu = (1-beta_I*(1+(1-tau_y)*r))/(1-beta_I*rho_b);

% <Capital and housing producers Block>
% (23) Investment in capital:
q_k = 1;
% (24) Law of motion of capital:
delta_k = ik_y/(k_y);
% (25) Investment in housing:
q_h = 1;
% (26) Law of motion of housing:
delta_h = ih_y/(h_y);

% (4) <Patient> Capital:
r_k = (1/beta_P - (1-delta_k) -tau_k*delta_k)/(1-tau_k);
% (3) <Patient> rental housing:
r_h = (1/beta_P - (1-delta_h-tau_p*(1-tau_y)) -tau_y*delta_h)/(1-tau_y);

% (20) Capital:
yn_y = r_k*(k_y)/alpha;
% (21) Utilization rate:
kappa_u = r_k;
% (35) Total investment:
i_y = ik_y + ih_y;
% (32) Goods market:
c_y = yn_y - i_y - g_y;
% (36) Definition of GDP:
rh_y = (1- (1+tau_c)*c_y- i_y - g_y);

% (12) Borrowing constraint:
hI_h = (b_h) * (1/phi) ;
% (33) Housing market:
hP_h = 1- hI_h - hR_h;
% (34) Total non-housing consumption:  (xi_h)
xi_h = ((h_y)/(c_y*(1+tau_c)))*(hP_h*(1-beta_P*(1-delta_h-tau_p*(1-tau_y))) + ...
    hI_h*(1-mu*(1-rho_b)*phi - beta_I*(1-delta_h-tau_p*(1-tau_y))) + hR_h*r_h);

% <Consumption block>
% (2) Patient C
cP_y = (hP_h)*(h_y)*(1-beta_P*(1-delta_h-tau_p*(1-tau_y)))/(xi_h*(1+tau_c));
% (8) Impatient C
cI_y = (hI_h)*(h_y)*(1-mu*(1-rho_b)*phi - beta_I*(1-delta_h-tau_p*(1-tau_y)))/(xi_h*(1+tau_c));
% (14) Renter C
cR_y = (hR_h)*(h_y)*r_h/(xi_h*(1+tau_c));
% Lambda*y
% (1) Patient Lambda
LambdaP = 1/((1+tau_c)*cP_y);
% (7) Impatient Lambda
LambdaI = 1/((1+tau_c)*cI_y);
% (13) Renter Lambda
LambdaR = 1/((1+tau_c)*cR_y);

% (Fiscal policy block) + (16; Renter BC):
Omega = tau_c*c_y + tau_y*(1-alpha)*yn_y + ...
    tau_y*(r_h-delta_h)*hR_h*(h_y) + tau_k*(r_k-delta_k)*(k_y) + tau_b*r*(b_h)*(h_y) ...
    -tau_y*r*(b_h)*(h_y) + tau_p*(1-tau_y)*(h_y) - g_y -(chi_P+chi_I+chi_R)*(yn_y);
BC_R_0 = ((1+tau_c)*(cR_y) + r_h*(hR_h)*(h_y) - chi_R*yn_y)/((1-tau_R)*(1-alpha)*yn_y);
% (27)
bg_y = (Omega + (1-alpha)*yn_y*(tau_R-tau_y)*BC_R_0)/((1-tau_b)*r-3*rrho_b - (tau_R-tau_y)*rrho_b/(1-tau_R));
% (31)
tax_y = (r-3*rrho_b)*bg_y + g_y +(chi_P+chi_I+chi_R)*(yn_y);
% (28-30)
trP_y = chi_P*yn_y - rrho_b*bg_y;
trI_y = chi_I*yn_y - rrho_b*bg_y;
trR_y = chi_R*yn_y - rrho_b*bg_y;

% (16): Renter BC
psi_R = BC_R_0 + (rrho_b/((1-alpha)*(1-tau_R)*yn_y))*bg_y;
% (11): Impatient BC
BC_I_LHS = (1+tau_c)*(cI_y) + (delta_h+tau_p*(1-tau_y))*(hI_h)*(h_y) + r*(1-tau_y)*(b_h)*(h_y);
psi_I = (BC_I_LHS - chi_I*yn_y + rrho_b*bg_y)/((1-alpha)*(1-tau_y)*yn_y);
% (CRTS)
psi_P = 1-psi_I-psi_R;

% <Labor Block>
% (5): Patient labor decision
xi_l = psi_P*(1-tau_y)*(1-alpha)*yn_y/((l_P)^(theta+1)*(1+tau_c)*cP_y);
% (9): Impatient labor decision
l_I = (psi_I*(1-tau_y)*(1-alpha)*yn_y/(xi_l*(1+tau_c)*cI_y))^(1/(theta+1));
% (15): Renter labor decision
l_R = (psi_R*(1-tau_R)*(1-alpha)*yn_y/(xi_l*(1+tau_c)*cR_y))^(1/(theta+1));
% (17): wage P
wP_y = (1-alpha)*psi_P*yn_y/l_P;
% (18): wage I
wI_y = (1-alpha)*psi_I*yn_y/l_I;
% (19): wage R
wR_y = (1-alpha)*psi_R*yn_y/l_R;

% Production
z_y = yn_y/((k_y)^alpha*(l_P^psi_P*l_I^psi_I*l_R^psi_R)^(1-alpha));



