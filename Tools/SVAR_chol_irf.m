function irfs=SVAR_chol_irf(Y, p, B, Sigma, horizon)
% Y라는 data를 받고, p라고 하는 시차를 입력, 
% B와 Sigma는 "reduced_VAR_OLS" 함수의 축약형 회귀 결과물
% horizon는 IRF 생성의 horizon
% 촐레스키 ordering에 의한 recursive로 SR SVAR을 계산하여
% 다음의 공식에 따라 (orthogonalized) IRF를 계산
% Psi_h = B1 * Psi_(h-1) + B2 * Psi_(h-2) + ... + Bp * Psi_(h-p) if h>p
% 에 따라 irfs라고 하는 충격반응함수 결과물을 return한다.
    [~, K]=size(Y);
    % 촐레스키 분해. matlab에서는 transpose를 취해주어야 함!
    A = chol(Sigma)';
    % 시차별로 회귀계수를 나누어서 설정하며, 여기서 end-1은 상수항을 버린다는 의미
    B_p = reshape(B(1:(end-1),:)',K,K,p);
    % irf값이 들어갈 빈칸 생성
    irfs = zeros(K,K,horizon);
    
    irfs(:,:,1) = A;    % impact period (h=0)
    for h=2:horizon     % h = 1,..,(horizon-1)
        for j=1:min(p,h-1)
            irfs(:,:,h) = irfs(:,:,h) + B_p(:,:,j)*irfs(:,:,h-j);
        end
    end
end