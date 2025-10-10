function [stat_dist] = stationary(P)
    tol = 1e-6;
    [V, D] = eig(P');
    lambda = diag(D);
    idx = find(abs(lambda - 1) < tol);
    eig_vector = V(:,idx);
    stat_dist = eig_vector/(sum(eig_vector));
    % 검산
    error = sum(abs(stat_dist'*P - stat_dist'));
    if error > tol
        fprintf("경고! stationary distribution이 항등식을 만족하지 않습니다.")
    end
end