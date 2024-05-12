% Core and Periphery(2024)
% Yeongwoong Do

% Step01
% Calculating the implied expected interest rate in 3Y KTB and
% Rescaleing the policy news shocks
clear
clc

%% Data
% 2008:8 ~ 2024:2
data_KTB = readmatrix('MP_time.xlsx','Sheet','data','Range','A3:F162');
data_KOSPI = readmatrix('MP_time.xlsx','Sheet','data','Range','J3:L162');
[T, K] = size(data_KTB);

%% (1) Calculation the implied interest rate
% parameter
C1 = 8;      % coupon rate ~2010.9
C2 = 5;      % coupon rate 2010.10~
V = 100;    % face value

% 2010.10  
change_time = 28;

date = data_KTB(:,1);
Rate=zeros(T,2);
options=optimset('Display','off','TolX',1e-10,'TolFun',1e-10);
for j=1:2
    for t=1:T
        F = data_KTB(t,4+j);
        if t < change_time
            future_price = @(R) F - (C1/2)/(1+R/2) - (C1/2)/((1+R/2)^2) - (C1/2)/((1+R/2)^3) - (C1/2)/((1+R/2)^4) -(C1/2)/((1+R/2)^5) - ((C1/2)+100)/((1+R/2)^6);
        else
            future_price = @(R) F - (C2/2)/(1+R/2) - (C2/2)/((1+R/2)^2) - (C2/2)/((1+R/2)^3) - (C2/2)/((1+R/2)^4) -(C2/2)/((1+R/2)^5) - ((C2/2)+100)/((1+R/2)^6);
        end
        Rate(t,j) = fsolve(future_price,0.035,options)*100;  % percent(%)
    end
end

delta_Rate = Rate(:,2)-Rate(:,1);
data_KTB = [data_KTB Rate delta_Rate];

% ploting the policy news shocks
%figure(1)
%clf;
%plot(delta_Rate,'Linewidth',2)

%% (2) Regression for rescaling
% Refer to the Table 1 in Nakamura & Steinsson(2018, QJE)

market_rate= readmatrix('MP_time.xlsx','Sheet','data','Range','O3:AE162');
[~, R_num] = size(market_rate);
X = [ones(T,1) delta_Rate]; % percent(%)
beta_vec = zeros(R_num,1);
se_beta_vec = zeros(R_num,1);

for i=1:R_num
    Y = market_rate(:,i);
    beta = X'*X \ X'*Y;
    e = Y - X*beta;
    se_beta = ((e'*e)/(T-2))*inv(X'*X);
    beta_vec(i) = beta(2);
    se_beta_vec(i) = se_beta(2,2);     
end

% rescaling target = KORIBOR 100bp = 1 news shock
rescale = (beta_vec(2));
beta_vec_re = beta_vec./rescale;
se_beta_vec_re = se_beta_vec./rescale;
rescale_surprise = delta_Rate*rescale;

% check
% Y = market_rate(:,2);
% X = [ones(T,1) rescale_surprise];
% beta = X'*X \ X'*Y

%% (3) Checking for the existence of the information channel
% Refer to the Figure 1 in Jarocinski and Karadi(2020, AEJ_mac)
S_rate = delta_Rate;
S_stock = data_KOSPI(:,3); 

figure(1);
clf;
scatter(S_rate, S_stock,"filled")
xline(0)
yline(0)
xlabel("Surprise in the three-year KTB futures")
ylabel("Surprise in the KOSPI 200 futures")

%% (4) Organize the surprise series as monthly data

% convert excel date to matlab date
date = datetime(date, 'ConvertFrom', 'excel');

% set the time
start_date = datetime('2008-08-01');
end_date = datetime('2024-02-01');
month_vector = (start_date:calmonths(1):end_date)';

% fill the data
month_surprise_rate = zeros(length(month_vector),1);
month_surprise_stock = zeros(length(month_vector),1);
for i = 1:T
    % extract year and month
    [year, month, day, ~, ~, ~] = datevec(date(i));
    % location in month_vector
    if year==2008
        idx = (month-7);
    elseif year > 2008
        idx = 5 + (year-2009)*12 + month;
    end
    
    % 해당 월의 데이터 누적
    month_surprise_rate(idx) = month_surprise_rate(idx) + S_rate(i);
    month_surprise_stock(idx) = month_surprise_stock(idx) + S_stock(i);
end

% save as csv
save('surprise.mat',"month_vector","month_surprise_rate", "month_surprise_stock");
writematrix([month_surprise_rate, month_surprise_stock], 'mp_shock_monthly.csv');
