function [beta, se0, se_r] = OLS_White(Y, X_nocons)
    T = length(Y);
    X = [X_nocons, ones(T,1)];
    XX = X'*X;
    beta = XX \ X'*Y;
    e = Y - X*beta;
    K = size(X, 2);
    
    % homoskedasticity
    V0 = inv(XX)*((e'*e)/(T-K));
    
    % White (heteroskedasticity-robust) variance-covariance matrix
    Omega = diag(e.^2);
    V_r = (T/(T-K))*inv(XX)* X'*Omega*X * inv(XX);

    se0 = sqrt(diag(V0));
    se_r = sqrt(diag(V_r));
end