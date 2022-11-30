function [pi_optimal, f_m_optimal, k_prime_c_optimal, alpha_0_optimal, alpha_1_optimal, I_0_optimal ] = LMW_0921(d, lambda, tau_m, tau_c, kappa_m, kappa_c)

%%% The static model replication of the  Li, McAndrew, Wang(2020,JME)
%%% d; unit cost of cardpayment
%%% lambda; 1/lambda = mean income
%%% tau_c; unit cost of trading cash for consumer
%%% tau_m; unit cost of trading cash for merchants
%%% k_m; fixed cost of trading card for consuper / mean income
%%% k_c: fixed cost of trading card for merchants / mean income
%%% 

% Korea card fee
k_prime_m = 0.00;
f_c = 0.00;

%grid
[f_m_grid,k_prime_c_grid] = ndgrid(d:0.0001:0.0206, -kappa_c:0.000001:0.0015);
[N, M]=size(k_prime_c_grid);


% profit function grid search
pi=zeros(N,M);
alpha_0=zeros(N,M);
alpha_1=zeros(N,M);
I_0=zeros(N,M);
for i=1:N
    for j=1:M
        k_m=k_prime_m+kappa_m;
        k_c=k_prime_c_grid(i,j)+kappa_c;
        f_m=f_m_grid(i,j);
        Z_1=((1-f_m)/(1+f_c)-(1-tau_m)/(1+f_c));
        Z_0=((1-f_m)/(1+f_c)-(1-tau_m)/(1+tau_c));
        % Boundary
        LB_alpha_0 = max( [ k_m/(Z_0*2*(1-k_c)), 0]);
        UB_alpha_0 = 1;
        LB_I_0 = max([ (k_c*(1+tau_c))/((1+f_c)*lambda) / ((1+tau_c)/(1+f_c)-1), 0]);
        % initial guess = LB_I_0
        I_old=LB_I_0;
        % prepare for loop
        tolerance=10^(-6); diff = 10;
        iter=0; max_iter=1000;
        while ((diff>tolerance) && (iter<=max_iter))
                alpha_next = k_m./(2*exp(-lambda.*I_old).*(1+lambda.*I_old-k_c).*Z_0);
                alpha_1_in = (Z_0/Z_1)*alpha_next;
                up = min( [alpha_1_in, 1],[],2);
                A_term = intpart(up, tau_m, f_m, f_c, alpha_next, Z_0)...
                        - intpart(alpha_next, tau_m, f_m, f_c, alpha_next, Z_0);
                I_next= (((1+tau_c)/(1+f_c))^(1-alpha_next^2)*(k_c/lambda))...
                    / (((1+tau_c)/(1+f_c))^(1-alpha_next^2) - exp(A_term));  
                diff = abs(real(I_next)-real(I_old));
                if alpha_next >= UB_alpha_0
                    break
                end
                I_old = I_next;
                iter = iter + 1;
            end
            if ((alpha_next>=UB_alpha_0) || (alpha_next<LB_alpha_0) || (I_old<0) || (isreal(I_old)==0))
                pi(i,j)=NaN;
                alpha_0(i,j)=NaN;
                alpha_1(i,j)=NaN;
                I_0(i,j)=NaN;
            else
                pi(i,j)=((1-alpha_next^2)./(1+f_c)).*(exp(-lambda*I_old)./lambda)*(1+lambda.*I_old-k_c).*(f_c+f_m-d)+...
                    + (k_prime_c_grid(i,j)/lambda)*exp(-lambda*I_old);
                alpha_0(i,j)=alpha_next;
                alpha_1(i,j)=min((Z_0/Z_1)*alpha_next,1);
                I_0(i,j)=I_old; 
            end
    end
end

% profit maximization
format shortG
[obj, index]=max(pi,[],'all','linear');
if obj>=0
    pi_optimal=pi(index);
    f_m_optimal=f_m_grid(index);
    k_prime_c_optimal=k_prime_c_grid(index);
    alpha_0_optimal=alpha_0(index);
    alpha_1_optimal=alpha_1(index);
    I_0_optimal=I_0(index);
else % exit case
    pi_optimal=NaN;
    f_m_optimal=NaN;
    k_prime_c_optimal=NaN;
    alpha_0_optimal=NaN;
    alpha_1_optimal=NaN;
    I_0_optimal=NaN;
end

end

