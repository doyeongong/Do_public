function y = intpart(x, tau_m, f_m, f_c, alpha_0, Z_0)

%%%  indefinite intergral calculation in I_0 calculation
a = (1-tau_m);
b = (1-f_m);
c = alpha_0*(1+f_c)*Z_0;
y = c*((x/b) + ((c*log(b*x-c))/b^2)) + (x^2)*log((a*x)/(b*x-c));
end
