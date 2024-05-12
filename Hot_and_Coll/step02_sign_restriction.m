% step03 부호제약을 통한 전통적 MP 충격 추출
% Kilian & Lutkepohl Ch.13을 참고

% 필요파일: PCA.m, reduced_VAR_OLS.m 파일이 있어야 함.
clear 
clc
% 
tic

%% (1) News shock
% period : 2008-08 ~ 2024-02

% News shock
load('surprise.mat');
S_rate = month_surprise_rate;
S_stock = month_surprise_stock;

% 변수갯수(k) 
k = 2;
T = length(S_rate);

%% VAR
Y = [S_rate-mean(S_rate), S_stock-mean(S_stock)]; 
Sigma = Y'*Y;

% 촐레스키 분해
A = chol(Sigma,'lower');

% P행렬 만들기 in (0, 2π)
N = 5000;
LB = 0;
UB = 2 * pi;
% (0, 2π) 범위에서 난수 생성
rng(20230812)
A_candidate = zeros(k,k,N);
e_candidate = zeros(T,k,N);
count = 1;

% P계산후 맞는지 검정
while count<N
    omega = LB + (UB - LB) * rand(1,1);
    P_candidate= [ cos(omega), -sin(omega); 
                   sin(omega), cos(omega)];
    A_tilde = A * P_candidate;
    if A_tilde(1,1)>0 && A_tilde(2,1)<0 && A_tilde(1,2) >0 && A_tilde(2,2)>0
        A_candidate(:,:,count) = A_tilde;
        for t=1:T
            e_candidate(:,:,count) = (A_tilde \ Y')';
        end
        count = count + 1;
    end
end
e_pc_seq = squeeze(e_candidate(:,1,1:count-1));
e_nc_seq = squeeze(e_candidate(:,2,1:count-1));

% median
e_tr = median(e_candidate(:,1,1:count-1), 3);
e_info = median(e_candidate(:,2,1:count-1), 3);
% mean
e_tr_mean = mean(e_candidate(:,1,1:count-1), 3);
e_info_mean = mean(e_candidate(:,2,1:count-1), 3);
% pca
[e_tr_pca, ~, ~, ~] = PCA(e_pc_seq, 1);
[e_info_pca, ~, ~, ~] = PCA(e_nc_seq, 1);

%% Meaning of the MP shocks scale

% Refer to the Table 1 in Nakamura & Steinsson(2018, QJE)

e_MP = e_tr_pca(1:102);
e_MP = [e_MP; e_tr_pca([103,105,106,108,109,111,112,114,115,117,118,120,121,123,124,126,127,129,130,132,133,135,136,138,139,141,142,144,145,...
    147,148,150,151,153,154,156,157,159,160,162,163,165,166,168,169,171,172,174,175,177,178,180,181,183,184,186,187])]; 

market_rate= readmatrix('MP_time.xlsx','Sheet','data','Range','O3:AE162');
market_rate(3,:) = market_rate(3,:)+market_rate(4,:);
market_rate(4,:)=[];
[T,~]=size(market_rate);
[~, R_num] = size(market_rate);
X = [ones(T,1) e_MP]; % percent(%)
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
e_tr_pca_scale = e_tr_pca*rescale;

save('sign.mat','e_tr','e_info','e_tr_mean','e_info_mean','e_tr_pca','e_info_pca','e_tr_pca_scale');
% check
% Y = market_rate(:,2);
% X = [ones(T,1) rescale_surprise];
% beta = X'*X \ X'*Y

%% graph

figure(1);
clf;
subplot(1,2,1);

plot(month_vector,e_tr,'-r','Linewidth',1)
hold on
plot(month_vector,e_tr_mean,'--b','Linewidth',1)
plot(month_vector,e_tr_pca,'-k','Linewidth',1)
xlabel('year')
title('MP shocks')

subplot(1,2,2);
plot(month_vector,e_info,'-r','Linewidth',1)
hold on
plot(month_vector,e_info_mean,'--b','Linewidth',1)
plot(month_vector,e_info_pca,'-k','Linewidth',1)
xlabel('year')
title('CB information shocks')

legend('median','mean','PCA')

toc





