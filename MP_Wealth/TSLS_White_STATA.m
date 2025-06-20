function [beta, se] = TSLS_White_STATA(Y, X, Z, W)
    T = length(Y);
    R = [X, W, ones(T,1)];       % full regressors
    Z_all = [Z, W, ones(T,1)];   % full instruments
    K = size(R,2);

    % First stage projection matrix
    PZ = Z_all * ((Z_all' * Z_all) \ Z_all');
    X_hat = PZ*R;

    % 2SLS estimator
    beta = (R' * PZ * R) \ (R' * PZ * Y);
    u = Y - R * beta;

    % Middle term (Stata-style robust part)
    Omega = diag(u.^2);
    A = R' * PZ * R;
    V1 = (T/(T-K))*inv(A) * (X_hat'*Omega*X_hat) * inv(A);   
    V2 = (T/(T-K))*(X_hat'*X_hat) \ (X_hat'*Omega*X_hat) / (X_hat'*X_hat);

    % Full robust VCOV 
    % usually se1>se2
    se1 = sqrt(diag(V1));  
    se2 = sqrt(diag(V2));  
    
    se = (se1+se2)/2;
end
