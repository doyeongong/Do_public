function [b_h, se_b_h] = LP_cum(Y, S, H, M,I)
% H = horizon
% M = lags of AR
% S = tge shock variable: (Nx1)
% Y = dependent variable(log): (1+M+N+H)x1
%                   I
%                 |--|
% |------------|-----|------------|---------|
%       (x)        M       N           H
S_total = S;
S = S(I+1:end);
N = length(S);
b_h = zeros(H+1);
se_b_h = zeros(H+1);

% S와 같은 시점의 Y벡터
start_idx   = M+1+1+I;
end_idx     = M+N+1+I;
Y_t = Y(start_idx:end_idx);

% AR lags
if M~=0
    Y_lag = zeros(N,M);
    for m = 1:M
        % ln Y_(t-m) - ln Y_(t-m-1)
        Y_lag(:,m) = Y(start_idx-m:end_idx-m)-Y(start_idx-m-1:end_idx-m-1);
    end
end

% Shock lags
if I~=0
    S_lag = zeros(N,I);
    for i = 1:I
        % S_(t-i)
        S_lag(:,i) = S_total(I-i+1:I-i+N);
    end
end

% generate stocked Y
Y_sur = [];    
for h=0:H
        Y_temp = Y(start_idx+h:end_idx+h)-Y(start_idx-1:end_idx-1);
        Y_sur = [Y_sur; Y_temp];
end
Y_prime = reshape(Y_sur,N,H+1)';
Y_sur=reshape(Y_prime, (H+1)*N,1);
    
% generate X
X_sur = [];
constant_sur = repmat(eye(H+1),N,1);
S_sur = [];
S_lag_sur = [];
Y_lag_sur = [];
    for t=1:N
        S_t = [S(t)*eye(H+1)];
        S_sur=[S_sur; S_t];
        Y_lag_t = [];
        S_lag_t = [];
        
        if M~=0
        for m=1:M
            Y_lag_t_m = Y_lag(t,m)*eye(H+1);
            Y_lag_t = [Y_lag_t, Y_lag_t_m];
        end
        Y_lag_sur = [Y_lag_sur; Y_lag_t];
        end
        
        if I~=0
            for i=1:I
            S_lag_t_m = S_lag(t,i)*eye(H+1);
            S_lag_t = [S_lag_t, S_lag_t_m];
        end
        S_lag_sur = [S_lag_sur; S_lag_t];
        end
        
    end
    if M~=0 && I~=0
        X_sur = [constant_sur, S_sur, Y_lag_sur, S_lag_sur];
    elseif M~=0 && I==0
        X_sur = [constant_sur, S_sur, Y_lag_sur];
    elseif M==0 && I~=0
        X_sur = [constant_sur, S_sur, S_lag_sur];
    else
        X_sur = [constant_sur, S_sur];
    end
    
% estimation
[beta_sur, Var_DKSCC]=LP_with_Driscoll_Kraay(Y_sur,X_sur,H+1,12);
beta_Var = Var_DKSCC(H+2:2*H+2,H+2:2*H+2);
b_h = beta_sur(H+2:2*H+2);
se_b_h = sqrt(diag(beta_Var));
    
end
