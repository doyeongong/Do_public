function [pol_a_idx] = to_idx(pol_a, a_grid)
    [N, K] = size(pol_a);
    for n=1:N
        for k=1:K
            a_p = pol_a(n,k);
            [~, idx] = min(abs(a_p - a_grid));
            pol_a_idx(n,k) = idx;
        end
    end
end