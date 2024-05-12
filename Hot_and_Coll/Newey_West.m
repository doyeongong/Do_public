function [beta, Var_beta] = Newey_West(Y, X, m)
% Calcuate Newey-West estimator given Data and the number of lags(m)
% This function is exaclty same with STAT command "newey"
[n, k] = size(X);
beta = (X'*X)\(X'*Y);
e = Y - X*beta;
if m==0  % White Heteroskadasticity Consistent estimator
    X_Omega0_X = (n/(n-k))*X'*diag(e.^2)*X;
    Var_beta=inv(X'*X)*X_Omega0_X*inv(X'*X);
else
    X_Omega0_X = (n/(n-k))*X'*diag(e.^2)*X;
    Addition = 0;
    for l=1:m
        sum=0;
        for t=l+1:n
            sum= sum + e(t)*e(t-l)*(X(t,:)'*X(t-l,:)+X(t-l,:)'*X(t,:));
        end
        Addition = Addition + (n/(n-k))*(1-l/(m+1))*sum;
    end
    X_Omegam_X = X_Omega0_X + Addition;
    Var_beta=inv(X'*X)*X_Omegam_X*inv(X'*X);
end

end