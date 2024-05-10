function [beta_sur, Var_DKSCC]=LP_with_Driscoll_Kraay(Y_sur,X_sur,H,L)
% Estimating LP with SUR and Driscoll-Kraay SE

% reshape the data as SUR
% -----------------------
% t=1  region=1  grdp pop
% t=1  region=2  grdp pop
% ....
% t=1  region=7  grdp pop
%-------------------------
% t=2  region=1  grdp pop
% t=2  region=2  grdp pop
% ....
% t=2  region=7  grdp pop
%-------------------------
% H : horizon of the Local projection
% Y_sur : (TxH)x1
% X_sur : (T)x(HxK)
% beta_sur: HxK
% u_sur : (TxH)x1
% L : the number of lags in Driscall_Kraay SE
% Note: there is no small sample adjustment (same with the "ase" option! in STATA)

% timeseries size
[TH, HK] = size(X_sur);
T = TH/H;
H = HK/H;

% estimate beta
beta_sur = (X_sur'*X_sur)\(X_sur'*Y_sur);

% SE
u_sur = Y_sur - X_sur*beta_sur;
V_ct = 0;
for t=1:T
    u_t=u_sur((t-1)*H + 1: t*H );
    X_t = X_sur((t-1)*H+1:t*H,:);
    V_ct = V_ct + (X_t)'*u_t*u_t'*(X_t);
end
V_ctl = 0;
for l=1:L
    for t=1:T-l
        u_t = u_sur((t+l-1)*H+1:(t+l)*H);
        u_t_l = u_sur((t-1)*H+1:t*H);
        X_t = X_sur((t+l-1)*H+1:(t+l)*H,:);
        X_t_l = X_sur((t-1)*H+1:t*H,:);
        V_ctl = V_ctl + (1-l/(L+1))*( (X_t)'*u_t*u_t_l'*(X_t_l) + (X_t_l)'*u_t_l*u_t'*(X_t) );
    end
end

Omega_DKSCC = V_ct + V_ctl;
Var_DKSCC = ((X_sur'*X_sur) \ Omega_DKSCC) * inv(X_sur'*X_sur);

end


