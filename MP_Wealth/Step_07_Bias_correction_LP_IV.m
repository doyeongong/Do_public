%% Data
% Step 07 Bias correction LP IV
% Following Herbst and Johannsen(2024, JE)
% Yeongwoong Do
clear
clc

% lag and horizon 
diff = 1;  % Difference of Gini
lag = 4;   % number of lag in control variables
H   = 12;  % horizon of IRF

% load the data
wkdir = 'C:\Users\82106\Documents\5. 2025\MP_Wealth_inequality\replication_package';
filename = fullfile(wkdir, 'Table', 'T5_dataset.xlsx');
data = readtable(filename);

% index
% column location
Gini_idx = 15;
R_idx    = 2;
MPS_idx  = 6;

% row location
s_q1  = 7;   % 1989Q3
s_q2  = s_q1+lag+diff;  % 1990Q4 (lag + diff)
e_q1  = 144;  % 2023Q4
e_q2  = e_q1-H; 



%% LP_IV
Y = data{s_q2:e_q1,Gini_idx}-data{s_q2-1:e_q1-1,Gini_idx};
X = data{s_q2:e_q1,R_idx};
Z = data{s_q2:e_q1,MPS_idx};
N = length(X);

W = [];
for i=1:lag
    new_col = [data{s_q2-i:e_q1-i,R_idx}, data{s_q2-i:e_q1-i,Gini_idx}-data{s_q2-i-1:e_q1-i-1,Gini_idx}];
    W = [W, new_col];
end

% First stage
% Note: Stock & Yogo (2005)'s F statistics based on the homoskedasticity assumption
[beta_f, se_f0, se_fr] = OLS_White(X, [Z W]);
F_first = (beta_f(1)/se_f0(1))^2;

% main regression
beta_h = zeros(H+1,1);
se_h = zeros(H+1,3);
N_h  = zeros(H+1,2);
for h=0:H
    % Y, X, Z, W
    Y_h = data{s_q2+h:e_q1,Gini_idx}-data{s_q2-1:e_q1-h-1,Gini_idx};
    X_in   = data{s_q2:e_q1-h,R_idx};
    Z_in   = data{s_q2:e_q1-h,MPS_idx};
    N_h(h+1,1) = length(Y_h);
    W_in = [];
    for i=1:lag
        new_col = [data{s_q2-i:e_q1-h-i,R_idx}, data{s_q2-i:e_q1-h-i,Gini_idx}-data{s_q2-i-1:e_q1-h-i-1,Gini_idx}];
        W_in = [W_in, new_col];
    end
    [beta, se] = TSLS_White(Y_h, [X_in W_in], [Z_in W_in]);
    [~, se_s]  = TSLS_White_STATA(Y_h, X_in, Z_in, W_in);
    [~, se_nw] = TSLS_HAC_STATA(Y_h, [X_in W_in], [Z_in W_in], h);
    beta_h(h+1,1) = beta(1);
    se_h(h+1,1) = se(1);    
    se_h(h+1,2) = se_s(1);    
    se_h(h+1,3) = se_nw(1);     
    
end

IRF=zeros(H+1,7);
IRF(:,1)= beta_h(:,1);           % main line
IRF(:,2)= beta_h(:,1)-se_h(:,2); % lower bound
IRF(:,3)= beta_h(:,1)+se_h(:,2); % upper bound
IRF(:,4)= beta_h(:,1)-1.645*se_h(:,2);
IRF(:,5)= beta_h(:,1)+1.645*se_h(:,2);
IRF(:,6)= beta_h(:,1)-se_h(:,3); 
IRF(:,7)= beta_h(:,1)+se_h(:,3); 

%% Bias correction

T = size(W,1);

% ACF term in bias correction
acf_corr = nan(1,H);
w = W - mean(W); % de-mean
Sigma_0 = cov(w); % var-cov
    for j=1:H
        Sigma_j = (w(1:end-j,:)'*w(j+1:end,:))/(T-j-1); % j-th autocov
        acf_corr(j) = 1+trace(Sigma_0\Sigma_j);
    end

% Iterate on bias correction
IRF_corr = IRF(:,1);
for h=1:H
    IRF_corr(h+1) = IRF(h+1) + (1/(T-h))*acf_corr(1:h)*IRF_corr(h:-1:1);
end


%% IRF figure

timeline=(0:1:H)';

figure(1)
clf

x_fill = [timeline; flipud(timeline)];
y_fill = [IRF(:,3); flipud(IRF(:,2))];
fill(x_fill, y_fill, [1 0 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

hold on
y_fill2 = [IRF(:,5); flipud(IRF(:,4))];
fill(x_fill, y_fill2, [1 0 0], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
plot(timeline,zeros(H+1,1),'--k', 'Linewidth', 1)  
plot(timeline,IRF(:,1),'-r', 'Linewidth', 2)
plot(timeline,IRF_corr,'--r', 'Linewidth', 2)

nw_lb = plot(timeline,IRF(:,6),':r', 'Linewidth', 2);
nw_lb.Color(4) = 0.3; 
nw_ub = plot(timeline,IRF(:,7),':r', 'Linewidth', 2);
nw_ub.Color(4) = 0.3; 

hold off
xlabel("Quarter")
ylabel("Difference of Gini")

figname = fullfile(wkdir, 'Figure', 'F3_main.png');
saveas(gcf, figname)   



figure(2)
clf

x_fill = [timeline; flipud(timeline)];
y_fill = [IRF(:,3); flipud(IRF(:,2))];
fill(x_fill, y_fill, [1 0 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

hold on
plot(timeline,zeros(H+1,1),'--k', 'Linewidth', 1)  
plot(timeline,IRF(:,1),'-r', 'Linewidth', 2)
nw_lb = plot(timeline,IRF(:,6),':r', 'Linewidth', 2);
nw_lb.Color(4) = 0.3; 
nw_ub = plot(timeline,IRF(:,7),':r', 'Linewidth', 2);
nw_ub.Color(4) = 0.3; 

hold off
xlabel("Quarter")
ylabel("Difference of Gini")

figname = fullfile(wkdir, 'Figure', 'F3_White_NW.png');
saveas(gcf, figname) 


figure(3)
clf

x_fill = [timeline; flipud(timeline)];
y_fill = [IRF(:,3); flipud(IRF(:,2))];
fill(x_fill, y_fill, [1 0 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

hold on
y_fill2 = [IRF(:,5); flipud(IRF(:,4))];
fill(x_fill, y_fill2, [1 0 0], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
plot(timeline,zeros(H+1,1),'--k', 'Linewidth', 1)  
plot(timeline,IRF(:,1),'-r', 'Linewidth', 2)
plot(timeline,IRF_corr,'--r', 'Linewidth', 2)

hold off
xlabel("Quarter")
ylabel("Difference of Gini")

figname = fullfile(wkdir, 'Figure', 'F3_Bias_corrected.png');
saveas(gcf, figname) 


% Table
IRF_all = [IRF, IRF_corr];  % Note: (H+1) x 8 
var_names = {'beta', 'w_lb', 'w_ub', 'w_lb_90', 'w_ub_90', 'nw_lb', 'nw_ub', 'bias_corrected'};
IRF_table = array2table(IRF_all, 'VariableNames', var_names);
IRF_table.Quarter = timeline; 
IRF_table = movevars(IRF_table, 'Quarter', 'Before', 1);

filename_out = fullfile(wkdir, 'Table', 'T3_main.xlsx');
writetable(IRF_table, filename_out);


