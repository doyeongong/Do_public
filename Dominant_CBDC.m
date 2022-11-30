function [pi_optimal, f_m_optimal, k_prime_c_optimal, alpha_1_optimal, alpha_0_optimal, I_0_optimal ] = Dominant_CBDC(d, lambda, y_m, y_c, kappa_m, kappa_c)

% rho >= tau_c 인 경우. 즉 CBDC가 dominant한 경우의 균형 계산

k_prime_m=0;
k_m = k_prime_m+kappa_m;
I_0_optimal=0;
k_c = 0;
k_prime_c_optimal = -kappa_c;
f_c = 0;

% FOC of the card network profit maximization
syms fm
sol=vpasolve(-2*k_m^2*(fm-d)+4*(y_m-fm)^3 -k_m^2*(y_m-fm)==0, fm, [0 y_m-0.000001]);
if sol<d
    f_m_optimal=d;
    pi_optimal=0;
else
    f_m_optimal=sol;
end

% alpha_1_optimal
alpha_1_optimal= min(round(k_m/(2*(1-k_c)*(y_m-f_m_optimal)),5),1);
Z_1=((1-f_m_optimal)/(1+f_c)-(1-y_m)/(1+f_c));
Z_0=((1-f_m_optimal)/(1+f_c)-(1-y_m)/(1+y_c));
if (Z_1/Z_0)>0
    alpha_0_optimal= min( (Z_1/Z_0)*alpha_1_optimal ,1);
else 
    alpha_0_optimal=1;
end

% profit
pi_optimal= max(round(((1-alpha_1_optimal^2)/(1+f_c)) * (1/lambda)*(1-k_c)*(f_m_optimal+f_c-d)+k_prime_c_optimal*(1/lambda),5),0);

end