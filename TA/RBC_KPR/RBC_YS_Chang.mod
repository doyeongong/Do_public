// RBC example 1 (without a labor choice)
// Yeongwoong Do, July 2024

// Endogeneous variables (10 variables)
var y k a c i n w r b lambda;
// Exogenous variables (2 shocks)
varexo zeta vareps;

// Parameters
parameters eis gamma alpha rho bar_r g; 
parameters gamma_x delta MPK eta I_Y;
eis   = 1;             % elasticity of intertemporal substitution
gamma = 2;           % Frisch elasticities
alpha = 0.66;        % labor income share
rho   = 0.75;        % autocorrelation coefficient for ln B
bar_r = 0.02;       % SS interest rate
g     = 0.00;       % ss net growth rate
gamma_x = exp(g);    % ss gross growth rate
delta = 0.05;      % depreciation
MPK   = bar_r + delta;
eta   = MPK/(MPK+(1-delta));
K_Y   = ((1-alpha)/MPK);
I_Y   = (gamma_x-(1-delta))*(K_Y);
						
// Model equations
model;
// (1) FOC for [c]
-(1/eis)*c = lambda;
// (2) FOC for [n]
b + gamma*n = lambda + w;
// (3) FOC for [k]
lambda + vareps =  lambda(+1) + eta*(a(+1) - alpha*k + alpha*n(+1));
// (4) w = MPL
w = a + (1-alpha)*k(-1) + (alpha-1)*n;
// (5) r = MPK - delta
r = (MPK/(MPK-delta))*(a - alpha*k(-1) + alpha*n);
// (6) production
y = a + (1-alpha)*k(-1) + alpha*n;
// (7) aggregate demand
y = (1-I_Y)*c + I_Y*i;
// (8) law of motion
i = (gamma_x/(gamma_x-(1-delta)))*(k+vareps) -(1-delta)/(gamma_x-(1-delta))*k(-1);
// (9) AR(1) of ln B
b = rho*b(-1) + zeta;
// (10) A is redundant
a = 0;
end;


//Exgeneous Shocks 
shocks;
var vareps = 0.01^2;
end;

model_diagnostics;
steady;

//Creating Impulse Response
stoch_simul(order = 1, irf=20, nocorr, nodecomposition, nofunctions, nomoments);