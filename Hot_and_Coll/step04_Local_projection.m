% Core and Periphery(2024)
% Yeongwoong Do

% Step04
% Local Projection
% cumulative LP (same with GKK(2023))
% ln Y_{j,t+h} - ln Y_{j,t-1} = alpha(h) + beta(h) S_t + AR terms + e_{j,t}
clear
clc
tic;
timeStart=datestr(clock);

%% Control

% region on/off
region=1;

% store results on/off
store=0;

% Horizon of IRFs (30개월 후까지, 31개의 점을 만듦)
H = 24; 
% M = lags of AR
M = 12;
% I = lags of shocks
I = 1;

% Total number of samples
% S_t 갯수 - H
% N = length(S)-H;
pre_del = 0;
% N = 129; % 2009.4~2019.12 (129)
N = 163; % 2008.8~2022.2 (163)
% N = 144; % 2008.8~2020.7 (144)

% MP shock을 무엇으로 할건지
% News shock
load('surprise.mat');
MP = month_surprise_rate;
SP = month_surprise_stock;

% sign restriction 
load('sign.mat');
load('OMPS.mat');
% S = e_tr_pca_scale;
% S = MP; 
S = e;

%% Loading Data

% period : 2008-08 ~ 2024-02: 187개

% macro variables : 2006-01 ~ 2024-02: 218개
macro = readmatrix('House_price_si_gun_gu.xlsx','Sheet','macro','Range','B9:J226');
[~, macro_num] = size(macro);
% sectoral house prices : 2006-01 ~ 2024-02: 218개
regional = readmatrix('House_price_si_gun_gu.xlsx','Sheet','real_trading_data','Range','B5:IM222');
[start_2006_1,J] =size(regional);

% deflate by CPI
for i=3:9
    if i~=4 && i~=5 && i~=6
        macro(:,i)=macro(:,i)./macro(:,2);
    end
end
for j=1:J
    regional(:,j) = regional(:,j)./macro(:,2);
end


 
%% Cumulative Local projection (Jorda(2005))

% 가져와야 하는 Y에서의 갯수
% M=0으로 AR시차변수가 없어도 과거 시차 한개가 필요함! (187+1개)
pre_period = start_2006_1 - pre_del - (M+N+H);
Y1 = log(macro(pre_period:end,:));
Y2 = log(regional(pre_period:end,:));
S = S(pre_del+1:N);

% estimation

fprintf('거시변수에 대한 Local projection을 시작합니다.\n')
b_mac = zeros(H+1,macro_num);
se_b_mac = zeros(H+1,macro_num);
for i=1:macro_num
    if i==5 || i==6
        [beta, std_beta] = LP_cum(macro(pre_period:end,i), S, H, M, I);
    else
        [beta, std_beta] = LP_cum(Y1(:,i), S, H, M, I);
    end
    b_mac(:,i)=beta;
    se_b_mac(:,i)=std_beta;
end

if region==1
fprintf('%d 개 지역별 주택가격에 대한 Local projection을 시작합니다.\n',J)
b_h = zeros(H+1,J);
se_b_h = zeros(H+1,J);
for j=1:J
    [beta, std_beta] = LP_cum(Y2(:,j), S, H, M, I);
    b_h(:,j)=beta;
    se_b_h(:,j)=std_beta;
    % 매 10번째 루프마다 실행 횟수 출력
    if mod(j, 10) == 0
        fprintf('현재까지 %d/246개 지역을 추정했습니다.\n', j);
    end
end
end

% store the results
if store==1
    results = struct;
    results.macro.data = macro;
    results.macro.beta = -b_mac; % to a expasionary MP shocks
    results.macro.se   = se_b_mac;
    results.y = regional;
    results.x = S;
    results.beta=-b_h; % to a expasionary MP shocks
    results.se = se_b_h;
    save('results.mat','results');
end

%% Ploting
 
fprintf('추정완료! IRF 그림그리기를 시작합니다.\n')
horizon = linspace(0,H,H+1);

% 거시변수
var_list = {'Output', 'CPI','Employment','Households debt','KORIBOR 3M','CD 91days','HP R-one','HP KB','HP REPS'};
figure(1);
clf;
sgtitle('IRF to expansionary MP shock')
for j=1:macro_num
     CI = [-b_mac(:,j)-1*se_b_mac(:,j), +2*se_b_mac(:,j)];
     subplot(3,3,j)
    area(horizon,CI,'FaceAlpha',0.3)
    newcolors = [1 1 1; 0 0.5 1];
    colororder(newcolors)
    hold on
     plot(horizon,-b_mac(:,j),'-r','Linewidth',2)
     plot(horizon,-b_mac(:,j)-2*se_b_mac(:,j),'--k','Linewidth',1)
     plot(horizon,-b_mac(:,j)+2*se_b_mac(:,j),'--k','Linewidth',1)
     plot(horizon,zeros(H+1,1),'-k','Linewidth',1)
    hold off
    title(var_list{j})
    if j==5 || j==6
        ylabel('percent(%)')
    else
        ylabel('log difference')
    end
    xlabel('month')
    xlim([0,H])
end

if region==1
% 지역별 IRF
figure(2);
clf;
subplot(1,2,1)
CI = [-b_mac(:,end)-1*se_b_mac(:,end), +2*se_b_mac(:,end)];
area(horizon,CI,'FaceAlpha',0.3)
newcolors = [1 1 1; 0 0.5 1];
colororder(newcolors)
hold on
plot(horizon,-b_mac(:,end),'-r','Linewidth',2)
plot(horizon,-b_mac(:,end)-2*se_b_mac(:,end),'--k','Linewidth',1)
plot(horizon,-b_mac(:,end)+2*se_b_mac(:,end),'--k','Linewidth',1)
plot(horizon,zeros(H+1,1),'-k','Linewidth',1)
hold off
ylim([-1 1])
xlim ([0, H])
ylabel('log difference')
xlabel('month')

subplot(1,2,2)
plot(horizon,zeros(H+1,1),'--k','Linewidth',1)
hold on
for j=1:J
    plot(horizon,-b_h(:,j),'-r','color',[0.5 0.8 1],'Linewidth',0.3)
end
plot(horizon,zeros(H+1,1),'--k','Linewidth',1)
plot(horizon,-b_mac(:,end),'-r','Linewidth',2)
hold off
ylim([-1 1])
xlim ([0, H])
ylabel('log difference')
xlabel('month')

end

%% Elapsed time

disp('--------------------------------------------------------')
fprintf('작업시작시각  : %s\n',timeStart)
fprintf('작업종료시각  : %s\n',datestr(clock))
toc;
disp('--------------------------------------------------------')




