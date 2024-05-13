% Core and Periphery(2024)
% Yeongwoong Do

% Analysis and Graph
clear
clc

% short name list of metropolitan areas
metro = {'NY', 'LA', 'CI','SD'};

% 루프를 사용하여 각 도시에 대한 데이터 로드 및 저장
for i = 1:length(metro)
    % 파일 경로 설정
    filePath = sprintf('C:\\Users\\82106\\Documents\\4. 2024학년도\\Core_Periphery\\2. Empirical\\US\\LP_result\\%s.mat', metro{i});
    % 데이터 로드
    data = load(filePath);
    % 결과를 구조체에 저장
    total.(metro{i}) = data.result;
end

%% Figure

% 2023m12
T = 288;
peak = 36;

long_name = {'New York', 'Los Angeles', 'Chicago','San Diego'};

for i=1:length(metro)
    price = log(total.(metro{i}).price(T,:))';
    Y = -total.(metro{i}).beta(peak,:)';
    X = [ones(length(price),1), price, price.^2];
    beta = (X'*X) \ (X'*Y);
    N_point = 100;
    price_dom = linspace(min(price),max(price),N_point)';
    x_dom = [ones(N_point,1), price_dom, price_dom.^2];
    Y_hat = x_dom*beta;

    subplot(2,2,i);
    scatter(price, Y,'filled')
    hold on
    plot(price_dom, Y_hat,'-r','Linewidth',2)
    hold off
    xlabel('log housing value')
    ylabel('IRF on 36 months later')
    title(long_name{i})
end

