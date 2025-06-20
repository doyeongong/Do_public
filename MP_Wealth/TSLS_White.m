function [beta, se] = TSLS_White(Y, X_nocons, Z_nocons)
    T = length(Y);
    X = [X_nocons, ones(T,1)];
    Z = [Z_nocons, ones(T,1)];
    X_hat = Z * ((Z' * Z) \ (Z' * X));
    beta = (X_hat'*X_hat) \ X_hat'*Y;
    e = Y - X*beta;
    K = size(X, 2);

    % White (heteroskedasticity-robust) variance-covariance matrix
    Omega = diag(e.^2);
    V = (T/(T-K))*(X_hat'*X_hat) \ (X_hat'*Omega*X_hat) / (X_hat'*X_hat);

    se = sqrt(diag(V));
end