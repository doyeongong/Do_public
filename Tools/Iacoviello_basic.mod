// The basic model of Iacoviello(2005, AER)
// Yeongwoong Do, June 2023

// endogeneous variables
var h1 Y c1 c2 b;
var R pi X_hat q rr;

// exogenous variables (shocks)
varexo e_R;

// Parameters
parameters m theta beta gamma nu j X eta r_y r_pi r_R;
parameters gamma_e n1 n2 bY qhY h_h iota kappa;

// Calibration
m 	= 0.89; 
theta 	= 0.75;
beta 	= 0.99;
gamma 	= 0.98;
nu 	= 0.03;
j 	= 0.1; 
X   = 1.05;
eta   = 1.01;
r_y    = 0;
r_pi 	= 0.27;
r_R = 0.73;

// calculation
gamma_e = (1-m)*gamma + m*beta;
n1 = nu*(1-gamma)*(1-beta*m)/((1-gamma_e)*X);
n2 = (X - nu + gamma*nu*(1-beta)*m/(1-gamma_e))*(1/X);
bY = beta*m*gamma*nu/((1-gamma_e)*X);
qhY = (gamma*nu)/((1-gamma_e)*X);
h_h = (gamma*nu*(1-beta))/(j*((X-nu)*(1-gamma_e)+gamma*nu*(1-beta)*m));
iota = (1-beta)*h_h;
kappa = (1-theta)*(1-beta*theta)/theta;


// Model equations
model(linear);
// (L1) total output
Y = n1*c1+ n2*c2;
// (L2) the Euler equation for household consumption
c2 = c2(+1) - rr;
// (L3) the entrepreneurial flow of funds
n1*c1 = bY*b+(1/beta)*bY*(pi - R(-1)-b(-1)) + (nu/X)*(Y-X_hat) - qhY*(h1 - h1(-1));
// (L4) consumption/housing margin for entrepreneurs
q = gamma_e*q(+1) + (1-gamma_e)*(Y(+1)-h1-X_hat(+1))-m*beta*rr - (1-m*beta)*(c1(+1)-c1);
// (L5) consumption/housing margin for HHs
q = beta*q(+1) +iota*h1 + c2 - beta*c2(+1);
// (L6) the borrowing constraint
b = q(+1) + h1 - rr;
// (L7) production function
Y = ((eta*nu)/(eta-(1-nu)))*h1(-1) - ((1-nu)/(eta-(1-nu)))*(X_hat + c2);
// (L8) NKPC
pi = beta*pi(+1) - kappa*X_hat;
// (L9) MP rule
R = (1-r_R)*((1+r_pi)*pi(-1)+r_y*Y(-1)) + r_R*R(-1) + e_R;
// (Def) rr
rr = R - pi(+1);
end;

//Exgeneous Shocks 
shocks;
var e_R = 0.29^2;
end;
   
//Check steady states
resid(1);
steady;
check;

//Creating Impulse Response
stoch_simul(order = 1,irf=40)
Y pi q R;


