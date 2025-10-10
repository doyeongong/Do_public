// RBC example 1 (without a labor choice)
// Yeongwoong Do, July 2024

// Endogeneous variables (7 variables)
var y k c n w R b;
// Exogenous variables (2 shocks)
varexo zeta;

// Parameters
parameters eis gamma alpha rho bar_r g; 
parameters gamma_x delta MPK eta I_Y K_Y;
eis   = 1;             % elasticity of intertemporal substitution
gamma = 2;           % Frisch elasticities
alpha = 0.66;        % labor income share
rho   = 0.75;        % autocorrelation coefficient for ln B
bar_r = 0.02;       % SS interest rate
g     = 0.00;       % ss net growth rate
gamma_x = exp(g);    % ss gross growth rate
delta = 0.05;      % depreciation
MPK   = bar_r + delta;
eta   = MPK/(MPK-(1-delta));
K_Y   = ((1-alpha)/MPK);
I_Y   = (gamma_x-(1-delta))*(K_Y);
						
// Model equations
model;
// (1) FOC for [n]
b + gamma*n = w - c;
// (2) FOC for [k]
-c  =  -c(+1) + R(+1);
// (3) w = MPL
w =  (1-alpha)*k + (alpha-1)*n;
// (4) r = MPK - delta
R = (MPK/(MPK+(1-delta)))*(- alpha*k + alpha*n);
// (5) production
y = (1-alpha)*k + alpha*n;
// (6) resource constraint
y = (1-K_Y+ K_Y*(1-delta))*c + K_Y*k(+1) - K_Y*(1-delta)*k;
// (7) AR(1) of ln B
b = rho*b(-1) + zeta;

end;


//Exgeneous Shocks 
shocks;
var zeta = 0.01^2;
end;

model_diagnostics;
steady;

//Creating Impulse Response
stoch_simul(order = 1, irf=20, nocorr, nodecomposition, nofunctions, nomoments);