function [beta, se] = TSLS_HAC_STATA(Y, X_nocons, Z_nocons, lag)
    % Two-stage least squares with Newey-West (HAC) standard errors
    % Y: dependent variable (T x 1)
    % X_nocons: endogenous regressors (T x K1)
    % Z_nocons: instruments (T x K2)
    % lag: Newey-West lag length (integer)

    T = length(Y);
    X = [X_nocons, ones(T,1)];
    Z = [Z_nocons, ones(T,1)];
    
    % First stage: project X on Z
    X_hat = Z * ((Z' * Z) \ (Z' * X));

    % Second stage: 2SLS beta
    beta = (X_hat' * X_hat) \ (X_hat' * Y);
    u = Y - X * beta;

    % Prepare for HAC (Newey-West) standard error
    K = size(X, 2);
    S = zeros(K, K);

    % score (moment conditions)
    scores = X_hat .* u;  % T x K matrix

    for h = 0:lag
        w_h = 1 - h / (lag + 1);  % Bartlett kernel
        Gamma_h = zeros(K, K);
        for t = h+1:T
            Gamma_h = Gamma_h + scores(t,:)' * scores(t - h, :);
        end
        if h == 0
            S = S + Gamma_h;
        else
            S = S + w_h * (Gamma_h + Gamma_h');  % add both lag and lead
        end
    end

    % HAC variance-covariance matrix
    XX_inv = inv(X_hat' * X_hat);
    V_HAC = (T / (T - K)) * XX_inv * S * XX_inv;
    V_HAC2 = (T / (T - K)) * (X_hat' * X_hat) \ S / (X_hat' * X_hat);

    se1 = sqrt(diag(V_HAC));
    se2 = sqrt(diag(V_HAC2));
    se  = (se1+se2)/2;
end
