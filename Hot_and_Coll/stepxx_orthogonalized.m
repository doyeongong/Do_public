% Core and Periphery(2024)
% Yeongwoong Do

% Orthogonalized MPS
clear
clc


%% Loading Data
% period : 2008-08 ~ 2024-02: 187개

% News shock
load('surprise.mat');
MP = month_surprise_rate;
T = length(MP);

% macro variables : 2006-01 ~ 2024-02: 218개
macro = readmatrix('House_price_si_gun_gu.xlsx','Sheet','macro','Range','B9:K226');
[~, macro_num] = size(macro);

% deflate by CPI
for i=3:9
    if i~=4 && i~=5 && i~=6
        macro(:,i)=macro(:,i)./macro(:,2);
    end
end
macro = [macro(:,1), macro(:,2), macro(:,4), macro(:,9), macro(:,10)];
% macro = [macro(:,1), macro(:,2), macro(:,9), macro(:,10)];
ln_macro = log(macro);

%% prepare for the regression

point_200808 = 32;
% lag variables
lag_num = 12;
[~, var_num] = size(macro);
X = ones(T, lag_num*var_num+1);
for var = 1:var_num
for lag = 1:lag_num
        % X(:,lag+(var-1)*lag_num+1) = ln_macro(point_200808-lag:end-lag,var); % log level
        X(:,lag+(var-1)*lag_num+1) = ln_macro(point_200808-lag:end-lag,var) - ln_macro(point_200808-lag-1:end-lag-1,var);
end
end

beta = (X'*X) \ (X'*MP);
MP_hat = X*beta;
e = MP - MP_hat;
R2 = 1- (e'*e)/(MP'*MP);

save('OMPS.mat','e');

%% Figure
% ploting
figure(2);
clf;
time = linspace(1,T,T);
plot(time, MP, '-b','Linewidth',2)
hold on
plot(time, e, '--r','Linewidth',2)
plot(time, zeros(T,1),'--k')
legend('MPS','Orthogonalized MPS','Location','best')

