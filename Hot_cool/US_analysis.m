% Hot and Cooling market
% Yeongwoong Do(2024)

% Step01
% Load the orthogonalized MPS suggested by Bauer and Swanson(2023)
clear
clc
tic;
timeStart=datetime("now");

disp('--------------------------------------------------------')
fprintf('작업시작시각  : %s\n',timeStart)

%% Data

% MP surprise: 1988m2~2019m12
OMPS = readmatrix('OMPS.csv');
OMPS = OMPS(:,2);

% Zillow datas: 2000m1~2024m3
HP_data = readtable('Zillow_zip_code.csv','ReadVariableNames', true);

% 'Metro' 매칭
rows = HP_data.Metro == "Denver-Aurora-Lakewood, CO";
% 해당하는 행에서 10번째부터 14번째 열까지의 데이터를 추출합니다.
data = table2array(HP_data(rows, 10:end))';
% NaN 포함한 열 확인
hasNaN = any(isnan(data), 1);
% NaN 포함 열 제거
data(:, hasNaN) = [];


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
Y = log(data(1:end-1,:));

% Horizon of IRFs 
H = 48; 
% M = lags of AR
M = 12;
% I = lags of shocks
I = 2;
% J = number of sectors
[~, J] = size(data);

b_h = zeros(H+1,J);
se_b_h = zeros(H+1,J);

for j=1:J
    [beta, std_beta] = LP_cum(Y(:,j), S, H, M, I);
    b_h(:,j)=beta;
    se_b_h(:,j)=std_beta;
    % 매 10번째 루프마다 실행 횟수 출력
    if mod(j, 10) == 0
        fprintf('현재까지 %d/%d개 지역을 추정했습니다.\n', j, J);
    end
end

result.beta = b_h;
result.std = se_b_h;
result.price = data;
save("DV.mat", "result");

%% Figure

fprintf('추정완료! IRF 그림그리기를 시작합니다.\n')

b_h = result.beta;
scatter(data(end-3,:),-b_h(37,:),"filled",'MarkerFaceColor', [0 0 1]);
ylabel('Impulse responses after 36 months')
xlabel('log house prices in December 2023')


fprintf('작업종료시각  : %s\n',datetime("now"))
toc;
disp('--------------------------------------------------------')
