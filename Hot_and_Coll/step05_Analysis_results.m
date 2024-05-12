% Core and Periphrey(2024)
% plotting
clear
clc

%% Load IRF results
load('results.mat');
region = results.y;
S = results.x;
% IRF under a expasionary MP shock
irfs = results.beta;
irfs_mac = results.macro.beta;
[irf_num, J] = size(irfs);
H = irf_num-1;
% Find peak horizon
[agg_peak, agg_peak_h] = max(irfs_mac(:,end));
% IRFS (peak of aggregate)
irfs_peak = irfs(agg_peak_h,:)';


%% Categorizing
% aggregate IRF, aggregate house prices
HP_2020 = readmatrix('실거래가_2020년평균.xlsx','Sheet','Sheet1','Range','H2:H249');
aggre_HP = HP_2020(end-1);
seoul_HP = HP_2020(end);
regional_HP = HP_2020(1:end-2);
med_irfs = median(irfs_peak);
cutoff = aggre_HP;

category = zeros(J,4);
for j=1:J
    % hot market (IRF above the aggregate IRF values)
    if irfs_peak(j)>=agg_peak
        category(j,1)=1;
    end
    % high cooling market
    if ~isnan(regional_HP(j)) && irfs_peak(j)<agg_peak && regional_HP(j)>cutoff
        category(j,2)=1;
    end
end
category(:,3)=ones(J,1)-category(:,1)-category(:,2);
category(:,4)=category(:,1)*2+category(:,2)*3+category(:,3)*1;

% generating a catgory index
hot_index = find(category(:,1)==1);
high_cooling_index = find(category(:,2)==1);
low_cooling_index = find(category(:,3)==1);

% categorical regional house price matrix
hot = region(:,hot_index);
high_cooling = region(:,high_cooling_index);
low_cooling = region(:,low_cooling_index);

%% PCA and get first principal component

% M = lags of AR
M = 12;
% Total number of samples
N = 144; % 2008.8~2020.7 (144)
% I = lags of shocks
I = 0;
pre_del = 0;
[start_2006_1,~] = size(region);

% 가져와야 하는 Y에서의 갯수
% M=0으로 AR시차변수가 없어도 과거 시차 한개가 필요함! (187+1개)
pre_period = start_2006_1 - (M+N+H);

% % Simple PCA
% [total, ~, ~, ~] = PCA(region, 1);
% [hot_pc, ~, ~, ~] = PCA(hot, 1);
% [high_cooling_pc, ~, ~, ~] = PCA(high_cooling, 1);
% [low_cooling_pc, ~, ~, ~] = PCA(low_cooling, 1);
% F=[hot_pc, high_cooling_pc, low_cooling_pc]+100;
% K=3;

% MLE-EM 
non_hot = find(category(:,1)~=1);
non_high = find(category(:,2)~=1);
non_low = find(category(:,3)~=1);
K=3;% 4개 (전체, hot, high_cool, low_cool)
[F, lambda, Psi] = MLE_EM_factor(region, K, non_hot, non_high, non_low);
F = F+100; % 음수에 로그 취함을 방지

% estimation
Y = log(F);
b_class = zeros(H+1,K);
se_b_class = zeros(H+1,K);
for i=1:K
    [beta, std_beta] = LP_cum(Y(:,i), S, H, M, I);
    b_class(:,i)=beta;
    se_b_class(:,i)=std_beta;
end


%% Figure

% A1. dot plot
figure(1);
clf;

subplot(2,3,1)
scatter(regional_HP,irfs_peak,"filled",'MarkerFaceColor', 'blue')
hold on
yline(0)
hold off
xlim([8.5,12.5])
ylabel('Impulse responses at the peak horizon')
xlabel('log house prices in 2023')


% A2. dot plot with quatratic regression line
% NaN이 포함된 행 제거
valid_indices = ~isnan(regional_HP);
% quadratic regression line
Y_valid = irfs_peak(valid_indices);
X_valid = [ones(length(Y_valid),1) regional_HP(valid_indices) regional_HP(valid_indices).^2];
beta = (X_valid'*X_valid)\(X_valid'*Y_valid);
x_dom = [ones(100,1), linspace(8.5,12.5,100)', linspace(8.5,12.5,100)'.^2];
y_hat = x_dom*beta;

subplot(2,3,2)
plot(x_dom(:,2), y_hat,'color',[0.5 0.8 1 0.3],'Linewidth',50)
hold on
scatter(regional_HP,irfs_peak,"filled",'MarkerFaceColor', 'blue')
yline(0)
hold off
xlim([8.5,12.5])
ylabel('Impulse responses at the peak horizon')
xlabel('log house prices in 2023')

% A3. SM vs others
SM = false(size(irfs_peak));
SM([1:69, 94:102]) = true;
SM_irf = irfs_peak(SM);
SM_HP = regional_HP(SM);
Others_irf = irfs_peak(~SM);
Others_HP = regional_HP(~SM);

subplot(2,3,3)
s =scatter(SM_HP,SM_irf,"filled",'MarkerFaceColor', [0 0 0]);
hold on
h = scatter(Others_HP,Others_irf,"filled",'MarkerFaceColor', [0 0 1]);
s.MarkerFaceAlpha = 0.2; % 마커 내부의 투명도 조절
h.MarkerFaceAlpha = 1; % 마커 내부의 투명도 조절
yline(0)
hold off
xlim([8.5,12.5])
ylabel('Impulse responses at the peak horizon')
xlabel('log house prices in 2023')
legend('Seoul metropolitan area','Other regions','location','best')

subplot(2,3,4)
s =scatter(SM_HP,SM_irf,"filled",'MarkerFaceColor', [0 0 1]);
hold on
h = scatter(Others_HP,Others_irf,"filled",'MarkerFaceColor', [0 0 0]);
s.MarkerFaceAlpha = 1; % 마커 내부의 투명도 조절
h.MarkerFaceAlpha = 0.2; % 마커 내부의 투명도 조절
yline(0)
hold off
xlim([8.5,12.5])
ylabel('Impulse responses at the peak horizon')
xlabel('log house prices in 2023')
legend('Seoul metropolitan area','Other regions','location','best')

subplot(2,3,5)
s =scatter(SM_HP,SM_irf,"filled",'MarkerFaceColor', [0 0 1]);
hold on
h = scatter(Others_HP,Others_irf,"filled",'MarkerFaceColor', [1 0 0]);
s.MarkerFaceAlpha = 0.5; % 마커 내부의 투명도 조절
h.MarkerFaceAlpha = 0.5; % 마커 내부의 투명도 조절
yline(0)
hold off
xlim([8.5,12.5])
ylabel('Impulse responses at the peak horizon')
xlabel('log house prices in 2023')
legend('Seoul metropolitan area','Other regions','location','best')



% A4. classification
subplot(2,3,6)
scatter(regional_HP(hot_index),irfs_peak(hot_index),"filled",'MarkerFaceColor', "red")
hold on
scatter(regional_HP(low_cooling_index),irfs_peak(low_cooling_index),"filled",'MarkerFaceColor', "blue")
scatter(regional_HP(high_cooling_index),irfs_peak(high_cooling_index),"filled",'MarkerFaceColor', "cyan")
yline(0)
hold off
xlim([8.5,12.5])
ylabel('Impulse responses at the peak horizon')
xlabel('log house prices in 2023')
legend('Hot market','Low cooling market','High cooling market','location','best')

% A5. heat map -> QGIS

% A6. PCA

horizon = linspace(0,H,H+1);
var_list = {'Hot market', 'High cooling market','Low cooling market'};
figure(2);
clf;
% sgtitle('IRF to expansionary MP shock')
for j=1:K
     CI = [-b_class(:,j)-1*se_b_class(:,j), +2*se_b_class(:,j)];
     subplot(1,3,j)
    area(horizon,CI,'FaceAlpha',0.3)
    newcolors = [1 1 1; 0 0.5 1];
    colororder(newcolors)
    hold on
     plot(horizon,-b_class(:,j),'-r','Linewidth',2)
     plot(horizon,-b_class(:,j)-2*se_b_class(:,j),'--k','Linewidth',1)
     plot(horizon,-b_class(:,j)+2*se_b_class(:,j),'--k','Linewidth',1)
     plot(horizon,zeros(H+1,1),'--k','Linewidth',1)
    hold off
    title(var_list{j})
    xlim([0, 24])
    ylim([-0.08 0.12])
    ylabel('log difference')
    xlabel('month')
end

% A7. simple average IRF

avg_hot = mean(irfs(:,hot_index)');
avg_high = mean(irfs(:,high_cooling_index)');
avg_low = mean(irfs(:,low_cooling_index)');
figure(3);
clf;
% sgtitle('IRF to expansionary MP shock')
    plot(horizon,avg_hot,'-r','Linewidth',2)
    hold on
     plot(horizon,avg_high,'-c','Linewidth',2)
     plot(horizon,avg_low,'-b','Linewidth',2)
     plot(horizon,zeros(H+1,1),'--k','Linewidth',1)
    hold off
    ylabel('percent(%)')
    xlabel('month')
    legend('Hot market','High cooling market','low cooling market', 'Location','best');




