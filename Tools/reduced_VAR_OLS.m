function [B, Sigma, U]=reduced_VAR_OLS(Y, p)
% Y라는 data를 받고, p라고 하는 시차를 입력하면
% 축약형 회귀에서의 B계수 행렬과
% Sigma라고 하는 축약형 오차 분산-공분산행렬을 연산
% Y_t = B_1 Y_{t-1} + B_2 Y_{t-2} + ... + B_p Y_{t-p} + B_0 
% B'=[B_1 B_2 .... B_p B_0] 형태
    [T, ~] = size(Y);
    % 종속변수
    y=Y(p+1:T,:);
    % 설명변수
    x=[];
    for i=1:p
        x = [x Y(p+1-i:T-i,:)];
    end
    X = [x ones(T-p,1)];
    B = (X'*X)\(X'*y); % 상수항이 맨 아래에 있는 beta 행렬
    U = y- X*B;
    Sigma = (U'*U)/(T-p);
end