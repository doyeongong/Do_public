% Hot and Cooling market
% Yeongwoong Do(2024)

% Step01
% Load the orthogonalized MPS suggested by Bauer and Swanson(2023)
clear
clc

%% Data

% MP surprise: 1988m2~2019m12
OMPS = readmatrix('OMPS.csv');
OMPS = OMPS(:,2);

% Zillow datas: 2000m1~2024m3
SD_HP = readmatrix('SD.xlsx','Sheet','SD_example','Range','B2:CR292');
% NaN 포함 여부 확인
hasNaN = any(isnan(SD_HP), 1);
% NaN 포함 열 제거
SD_HP(:, hasNaN) = [];


%% Local Projection

% H = horizon
% M = lags of AR
% S = tge shock variable: (N+I)x1
% Y = dependent variable(log): (M+N+H+1)x1
%                   I
%                 |--|
% |------------|-----|------------|---------|
%       (x)        M       N           H

% 2000m12~2019m12
S = OMPS(155:end);
% 2000m1~2024m2
Y = log(SD_HP(1:end-1,:));


% Horizon of IRFs 
H = 48; 
% M = lags of AR
M = 12;
% I = lags of shocks
I = 2;
% J = number of sectors
[~, J] = size(SD_HP);


b_h = zeros(H+1,J);
se_b_h = zeros(H+1,J);

for j=1:J
    [beta, std_beta] = LP_cum(Y(:,j), S, H, M, I);
    b_h(:,j)=beta;
    se_b_h(:,j)=std_beta;
end

result.beta = b_h;
result.std = se_b_h;

save("SD_results.mat", "result");

%%

result=load("SD_results.mat");
b_h = result.result.beta;
scatter(SD_HP(end-3,:),-b_h(37,:),"filled",'MarkerFaceColor', [0 0 1]);
ylabel('Impulse responses after 36 months')
xlabel('log house prices in December 2023')
