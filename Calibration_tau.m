function [tau_m_matched, tau_c_matched, min_distance ] = Calibration_tau(d, kappa_m, kappa_c)

%%% Korea data
%%% 2021년 지급수단 및 모바일금융서비스 이용행태 조사결과
%%% 카드소지자 비중(84.5%, 2021년 신용카드 기준)
%%% 카드거래비중 (74.4%, 2021년 금액기준, 현금, 계좌이체, 기타를 제외한 비중)
%%% f_m=2.06%로 고정시킨재 calibration 실행

%parameter
lambda=1/61.25;
f_m=0.0206; f_c=0.00;

%grid
[tau_m_grid,tau_c_grid] = ndgrid(0.025:0.0001:0.060,0.00:0.0001:0.030);
[N, M]=size(tau_m_grid);

% optimal value coresponding K gird
pi_vals=zeros(N,M);
k_prime_m_vals=zeros(N,M);
k_prime_c_vals=zeros(N,M);
alpha_0_vals=zeros(N,M);
I_0_vals=zeros(N,M);

% Loop 
for i=1:N
    for j=1:M
        tau_c=tau_c_grid(i,j);
        tau_m=tau_m_grid(i,j);
        [pi, k_prime_m, k_prime_c, alpha_0, I_0] = LMW_0825(d, lambda, tau_m, tau_c, kappa_m, kappa_c);
        pi_vals(i,j)=pi;
        k_prime_m_vals(i,j)=k_prime_m;
        k_prime_c_vals(i,j)=k_prime_c;
        alpha_0_vals(i,j)=alpha_0;
        I_0_vals(i,j)=I_0;
    end
end
        


Z_0_vals=((1-f_m)/(1+f_c)-(1-tau_m_grid)./(1+tau_c_grid));
Z_1_vals=((1-f_m)/(1+f_c)-(1-tau_m_grid)./(1+f_c));
alpha_1_vals=alpha_0_vals.*(Z_0_vals./Z_1_vals);

% target matching
target1=exp(-lambda.*I_0_vals); % share of card holders
card_trade=(1+lambda.*I_0_vals-k_prime_c_vals-kappa_c).*((1-(alpha_0_vals).^2)./(lambda.*(1+f_c))).*exp(-lambda.*I_0_vals);
cash_trade1=(1./(1+tau_c_grid)).*(1/lambda - (exp(-lambda.*I_0_vals)./lambda).*(1+lambda.*I_0_vals));
cash_trade2=(((alpha_0_vals).^2)./(1+tau_c_grid)).*((1+lambda.*I_0_vals-k_prime_c_vals-kappa_c).*(exp(-lambda.*I_0_vals))./lambda);
target2=card_trade./(card_trade+cash_trade1+cash_trade2); % share of card trading

distance=sqrt((target1-0.845).^2+(target2-0.836).^2);
[min_distance, index2]=min(distance,[],'all','linear');
tau_m_matched=tau_m_grid(index2);
tau_c_matched=tau_c_grid(index2);

end