% CBDC remuneration(rate)
clear

%parameter
d=0.005; lambda=1/61.25;
f_c=0.000; k_prime_m=0.00;
tau_m=0.0302; tau_c=0.0059; 
kappa_m=0.015; kappa_c=0.0005;

% tau experiment
CBDC_grid=0.0:0.001:0.020;
Q=length(CBDC_grid);
theta=0.0000;

pi_vec=zeros(Q,1);
f_m_vec=zeros(Q,1);
k_prime_c_vec=zeros(Q,1);
alpha_0_vec=zeros(Q,1);
alpha_1_vec=zeros(Q,1);
I_0_vec=zeros(Q,1);

for q=1:Q
    rho=CBDC_grid(q)+theta;
    y_c=tau_c-rho;  % y_c=tau_c;
    y_m=tau_m-rho;  % y_m=tau_m;

    if rho<(tau_c-f_c)
        [pi, f_m, k_prime_c, alpha_0, alpha_1, I_0] = LMW_0921(d, lambda, y_m, y_c, kappa_m, kappa_c);
        pi_vec(q,1)=pi;
        f_m_vec(q,1)=f_m;
        k_prime_c_vec(q,1)=k_prime_c;
        alpha_0_vec(q,1)=alpha_0;
        alpha_1_vec(q,1)=alpha_1;
        I_0_vec(q,1)=I_0;
    else
        [pi, f_m, k_prime_c, alpha_1, alpha_0, I_0] = Dominant_CBDC(d, lambda, y_m, y_c, kappa_m, kappa_c);
        pi_vec(q,1)=pi;
        f_m_vec(q,1)=f_m;
        k_prime_c_vec(q,1)=k_prime_c;
        alpha_0_vec(q,1)=alpha_0;
        alpha_1_vec(q,1)=alpha_1;
        I_0_vec(q,1)=I_0;
    end
end

% exit
alpha_0_vec(pi_vec==0)=NaN;
alpha_1_vec(pi_vec==0)=NaN;
I_0_vec(pi_vec==0)=NaN;
f_m_vec(pi_vec==0)=NaN;
k_prime_c_vec(pi_vec==0)=NaN;
k_c_vec=k_prime_c+kappa_c;

% Calculated index 
y_c_grid=tau_c-(CBDC_grid)'-theta*ones(Q,1);  % y_c_grid=tau_c-(rho_grid)';
y_m_grid=tau_m-(CBDC_grid)'-theta*ones(Q,1);  % y_m_grid=tau_m;
z1=((1-f_m_vec)./(1+f_c)-(1-y_m_grid)./(1+f_c));
z0=((1-f_m_vec)./(1+f_c)-(1-y_m_grid)./(1+y_c_grid));
alpha_threshold=min([alpha_0_vec, alpha_1_vec], [], 2);
s=log((1+f_c)./(1+y_c_grid))+(z0.*(1+f_c)./(1-f_m_vec)).*(z0./z1-1) ...
    +(z0.^2.*(1+f_c).^2./(1-f_m_vec).^2).*log((z0./z1-z0.*(1+f_c)./(1-f_m_vec))./(1-z0.*(1+f_c)./(1-f_m_vec)));

% share of card holder
share_card_holder=exp(-lambda.*I_0_vec);

% share of card trading
card_trade=(1+lambda.*I_0_vec-k_c_vec).*((1-(alpha_threshold).^2)./(lambda.*(1+f_c))).*exp(-lambda.*I_0_vec);
CBDC_trade1=(1./(1+y_c_grid)).*(1/lambda - (exp(-lambda.*I_0_vec)./lambda).*(1+lambda.*I_0_vec));
CBDC_trade2=(((alpha_threshold).^2)./(1+y_c_grid)).*((1+lambda.*I_0_vec-k_c_vec).*(exp(-lambda.*I_0_vec))./lambda);
CBDC_trade=CBDC_trade1+CBDC_trade2;
share_card_trade=card_trade./(card_trade+CBDC_trade);
share_card_trade(pi_vec==0)=0;
total_CBDC=(1/(1+y_c_grid))*(1/lambda);
CBDC_trade(pi_vec==0)=total_CBDC(pi_vec==0);
expense=CBDC_trade.*CBDC_grid'*2;

% consumer surplus
upper_I=exp(-lambda.*I_0_vec).*(lambda.*I_0_vec+1-k_prime_c_vec);
lower_I=1-exp(-lambda.*I_0_vec).*(1+lambda.*I_0_vec);
eta_1 = ((1+tau_c)*upper_I+(1+f_c)*lower_I)*(1-tau_m);
eta_2 = (1+tau_c).*(1-f_m_vec).*upper_I+(1+f_c)*(1-tau_m).*lower_I;
eta_3 = k_prime_m*0.5*(1+f_c)*(1+tau_c);

Phi_1=zeros(Q,1);
Phi_2=zeros(Q,1);
for q=1:Q
   e1=eta_1(q,1);
   e2=eta_2(q,1);
   e3=eta_3;
   Phi_1(q,1)=intpart(1, 1-e1, 1-e2, 1, 1, e3)-intpart(alpha_1_vec(q,1), 1-e1, 1-e2, 1, 1, e3);
   Phi_2(q,1)=intpart(alpha_1_vec(q,1), 1-e1, 1-e2, 1, 1, e3)-intpart(alpha_threshold(q,1), 1-e1, 1-e2, 1, 1, e3);
end
Phi_3=Phi_1+Phi_2;

Notcard=((1-y_m_grid)./(1+y_c_grid))*2.*exp(-Phi_1-0.5)*...
    (1/lambda).*(1-exp(-lambda.*I_0_vec)-lambda.*I_0_vec.*exp(-lambda.*I_0_vec));
card = ((1-y_m_grid)./(1+y_c_grid))*2.*exp(-Phi_3-0.5).*...
    ((1-y_c_grid)./(1+f_c)).^(1-alpha_threshold.^2).*(exp(-lambda.*I_0_vec)/lambda).*(1-k_prime_c_vec-kappa_c+lambda*I_0_vec);
CS=Notcard+card;

SW=CS-expense+pi_vec;

% exit case
CS_exit = (1/lambda).*((2*(1-y_m_grid))./(1+y_c_grid)).*exp(-0.5);
CS(pi_vec==0) = CS_exit(pi_vec==0);
expense_exit = (1/lambda).*CBDC_grid*2;
SW(pi_vec==0) = CS_exit(pi_vec==0) -expense_exit(pi_vec==0)';

% optimal rho
[SW_star, arg1]=max(SW,[],'all','linear');
rate_SW_star=CBDC_grid(arg1);
exit_index=(pi_vec==0);
exit_rate=CBDC_grid(find(exit_index,1));


output=[pi_vec, CS, SW, alpha_0_vec, alpha_1_vec, I_0_vec, f_m_vec, k_prime_c_vec, share_card_holder, share_card_trade];
xlswrite('output.xlsx',output);


%%
%%%%% graph

subplot(3,2,1);
plot(CBDC_grid,pi_vec,'linewidth',2)
title('Profit')
xlim([0,0.02])
xline(tau_c,'--')
xline(exit_rate,'--')
xlabel('\rho')
ylabel('\pi')


subplot(3,2,2);
plot(CBDC_grid,CS,CBDC_grid,SW,'linewidth',2)
xline(tau_c,'--')
xline(exit_rate,'--')
legend('CS','SW','','')
xlim([0,0.02])
title('Welfare')
xlabel('\rho')
ylabel('indirect utility')

subplot(3,2,3);
plot(CBDC_grid,alpha_0_vec,CBDC_grid,alpha_1_vec,'linewidth',2)
title('Thresholds for merchants')
xlim([0,0.02])
xline(tau_c,'--','color','k')
xline(exit_rate,'--','color','k')
xlabel('\rho')
ylabel('share')
ylim([0,1])
legend('\alpha_0','\alpha_1','Location','best')

subplot(3,2,4);
plot(CBDC_grid,I_0_vec,'linewidth',2)
title('Thresholds for consumers')
xlim([0,0.02])
xline(tau_c,'--','color','k')
xline(exit_rate,'--','color','k')
xlabel('\rho')
ylabel('Annual income')
ylim([0,40])

subplot(3,2,5);

yyaxis left
plot(CBDC_grid,f_m_vec,'color',[0.8500 0.3250 0.0980],'linewidth',2)
ylabel('rate(%)')

yyaxis right
plot(CBDC_grid,k_prime_c_vec,'color',[0 0.4470 0.7410],'linewidth',2)
ylabel('per annual income(%)')

title('Card fee')
xlim([0,0.02])
xline(tau_c,'--','color','k')
xline(exit_rate,'--','color','k')
xlabel('\rho')
ylabel('rate(%)')
legend('f_m','k_c','Location','best')

subplot(3,2,6);
plot(CBDC_grid,share_card_trade,CBDC_grid,share_card_holder,'linewidth',2)
title('Share of card')
xlim([0,0.02])
xline(tau_c,'--','color','k')
xline(exit_rate,'--','color','k')
xlabel('\rho')
ylabel('share(%)')
legend('card trading','card holder','Location','best')




