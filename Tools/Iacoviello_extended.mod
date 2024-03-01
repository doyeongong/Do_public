// The full model of Iacoviello(2005, AER)
// Yeongwoong Do, March 2024

// endogeneous variables
// 9 (same with the basic model)
var Y c c1 b h q pi R X;
// 8 (additional variables)
var c2 b2 h2 I K j u A;

// exogenous variables (4shocks)
varexo e_R e_j e_u e_A;

// Parameters
parameters beta gamma beta2 j_ss eta;
parameters mu nu psi delta phi_e phi_h X_ss theta;
parameters r_y r_pi r_R alpha m m2 rho_u rho_j rho_A sigma_u sigma_j sigma_A;
parameters zeta1 zeta2 zeta3 zeta4 gamma_e gamma_h;
parameters s1 s2 kappa b_ratio b2_ratio I_ratio c_ratio c1_ratio c2_ratio h_h1 h2_h1;

// Calibration
// 1. Preference 
beta 	= 0.99; 
gamma 	= 0.98;
beta2 	= 0.95;
j_ss 	= 0.1;
eta     = 1.01;
// 2. Technology
mu      = 0.3;
nu      = 0.03;
psi     = 2;
delta   = 0.03;
phi_e   = 0;
phi_h   = 0;
X_ss    = 1.05;
theta   = 0.75;
// 3. Policy stance
r_y     = 0.13;
r_pi 	= 0.27;
r_R     = 0.73;
// 4. Estimated
alpha   = 0.64;
m       = 0.89;
m2      = 0.55;
rho_u   = 0.59;
rho_j   = 0.85;
rho_A   = 0.03;
sigma_u = 0.17;
sigma_j = 24.89;
sigma_A = 2.24;

// calculation
zeta1 = gamma*mu / ((1-gamma*(1-delta))*X_ss);
gamma_e = (1-m)*gamma + m*beta;
zeta2 = gamma*nu / ((1-gamma_e)*X_ss);
zeta3 = j_ss/(1-beta);
gamma_h = beta2 + m2*(beta-beta2);
zeta4 = j_ss/(1-gamma_h);
b_ratio = m*beta*zeta2;
I_ratio = delta*zeta1;
c_ratio = (1- 1/beta)*b_ratio - I_ratio + (mu+nu)/X_ss;
s1 = (alpha*(1-mu-nu)+X_ss-1)/X_ss;
s2 = (1-alpha)*(1-mu-nu)/X_ss;
c2_ratio = s2/(1+(1-beta)*m2*zeta4);
c1_ratio = s1+ (1-beta)*(m*zeta2 + m2*zeta4*c2_ratio);
h_h1 = zeta2 / (zeta3*c1_ratio);
h2_h1 = (zeta4*c2_ratio) / (zeta3*c1_ratio);
b2_ratio = m2*beta*zeta4*(c2_ratio);
kappa = (1-theta)*(1-beta*theta)/theta;


// Model equations
model(linear);

// 1. Aggregate demand
// (A1) good market clearing condition
Y = c_ratio*c+ c1_ratio*c1 + c2_ratio*c2 + I_ratio*I;
// (A2) consumption Euler equation of patient HH
c1 = c1(+1) - R + pi(+1);
// (A3) optimal investment
I - K(-1) = gamma*delta*(I(+1)-K) + ((1-gamma*(1-delta))/psi)*(Y(+1)-X(+1)-K) + (1/psi)*(c-c(+1));

// 2. Housing/consumption margin
// (A4) entrepreneur
q = gamma_e*q(+1) + (1-gamma_e)*(Y(+1)-X(+1)-h)-m*beta*(R - pi(+1)) -(1-m*beta)*(c(+1)-c)- phi_e*(h-h(-1)-gamma*(h(+1)-h));
// (A5) impatient
q = gamma_h*q(+1) + (1-gamma_h)*(j-h2) - m2*beta*(R - pi(+1)) + (1-m2*beta)*c2 - (beta2-m2*beta2)*c2(+1) - phi_h*(h2-h2(-1)-beta2*(h2(+1)-h2));
// (A6) patient
q = beta*q(+1) + (1-beta)*j + (1-beta)*h_h1*h + (1-beta)*(h2_h1)*h2 + c1 - beta*c1(+1) + phi_h*(h_h1*(h-h(-1)) + h2_h1*h2 - beta*h_h1*(h(+1)-h) - beta* h2_h1*(h2(+1)-h2));

// 3. Borrowing constraint
// (A7) entrepreneur
b = q + h - R + pi(+1);
// (A8) impatient
b2 = q + h2 - R + pi(+1);

// 4. Aggregate supply
// (A9) production + labor market
Y = (eta/(eta-(1-mu-nu)))*(A + nu*h(-1) + mu*K(-1)) - ((1-mu-nu)/(eta-(1-mu-nu)))*(X + alpha*c1 +(1-alpha)*c2);
// (A10) PC curve
pi = beta*pi(+1) - kappa*X + u;

// 5. Flows of funds/evlolution of state variables
// (A11) lowa of motion for capital 
K = delta*I + (1-delta)*K(-1);
// (A12) entrepreneur's BC
b_ratio*b = c_ratio*c + zeta2*(h-h(-1)) + I_ratio*I + (b_ratio/beta)*(R(-1)+b(-1)-pi)-(1-s1-s2)*(Y-X);
// (A13) impatient's BC
b2_ratio*b2 = c2_ratio*c2 + (zeta4*c2_ratio)*(h2-h2(-1)) + (b2_ratio/beta)*(R(-1)+b2(-1)-pi)-s2*(Y-X);

// 6. Monetary policy and shock processes
// (A14) MP rule
R = (1-r_R)*((1+r_pi)*pi(-1)+r_y*Y(-1)) + r_R*R(-1) + e_R;
// (S1) preference shock
j = rho_j*j(-1) + e_j;
// (S2) markup shock
u = rho_u*u(-1) + e_u;
// (S1) TFP shock
A = rho_A*A(-1) + e_A;
end;

//Exgeneous Shocks 
shocks;
var e_R = 0.25^2;
end;
   
//Check steady states
resid(1);
steady;
check;

//Creating Impulse Response
stoch_simul(order = 1,irf=20)
Y pi q R;


