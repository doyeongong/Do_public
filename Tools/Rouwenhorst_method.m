function [Theta, y] = Rouwenhorst_method(rho, sigma, N)
    % print error message when N < 2
    if N < 2
        fprintf('Change N!\n');
    end
    
    p = (1 + rho) / 2;
    sigma_z = sigma / sqrt(1 - rho^2);
    psi = sqrt(N - 1) * sigma_z;
    
    % state space
    y_1 = -psi;
    y_N = psi;
    y = linspace(y_1, y_N, N);
    
    % transition matrix
    Theta = [p 1-p; 1-p p];
    if N == 2
        return;
    else
        Theta_old = Theta;
        for i = 3:N
            % Creating new transition matrix
            Theta_new = p * [Theta_old zeros(i-1, 1);  zeros(1, i)] + ...
                        (1-p) * [zeros(i-1, 1) Theta_old; zeros(1, i)] + ...
                        (1-p) * [zeros(1, i); Theta_old zeros(i-1, 1)] + ...
                        p * [zeros(1, i);  zeros(i-1, 1) Theta_old];
            row_sums = sum(Theta_new, 2); % normalizing by row sums
            Theta_old = Theta_new ./ row_sums;
        end
        Theta = Theta_old; % Update Theta to the latest Theta_old
    end
end