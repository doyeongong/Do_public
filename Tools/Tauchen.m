function [y_vec, Pi] = Tauchen(y_min, y_max, N, rho, sigma)
% Tauchen(1986, EL)에 의해 AR(1) 프로세스를 Markov transition matrix로 변형
% 논문의 y의 자기상관계수 lambda가 rho로 바뀌었음.

% y vector의 점 사이의 간의 간격
w=(y_max-y_min)/(N-1);
y_vec=zeros(1,N);
y_vec(1)=y_min;

% y vector 생성
for n=2:N
    y_vec(n)=y_vec(n-1)+w;
end

% Transition matrix Pi 생성(행:어제 -> 열:오늘)
Pi = zeros(N,N);
for j=1:N
    for k=1:N
        Left = (y_vec(k)-rho*y_vec(j)-w/2)/sigma;
        Right = (y_vec(k)-rho*y_vec(j)+w/2)/sigma;
        if k==1
            Pi(j,k)=normcdf(Right);
        elseif k==N
            Pi(j,k)=1-normcdf(Left);
        else
            Pi(j,k)=normcdf(Right)-normcdf(Left);
        end
    end
end

end