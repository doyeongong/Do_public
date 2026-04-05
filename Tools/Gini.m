function G = Gini(x, pop)
% compute_gini
% MATLAB version matched to the provided Fortran subroutine:
%   subroutine compute_gini(x, pop, N, G)
%
% Inputs:
%   x   : income/asset vector, length N
%   pop : population weight vector, length N
%
% Output:
%   G   : Gini coefficient

    % column vectors
    x = x(:);
    pop = pop(:);

    % size check
    N = length(x);
    if length(pop) ~= N
        error('x and pop must have the same length.');
    end

    % -----------------------
    % Step 1: argsort(x) -> idx
    % -----------------------
    [~, idx] = sort(x, 'ascend');

    % -------------------------------
    % Step 2: sorted xp(i) = x(i)*pop(i)
    % -------------------------------
    xp = x(idx) .* pop(idx);

    % -------------------------------
    % Step 3: total population and mean income
    % -------------------------------
    total_pop = sum(pop);
    total_income = sum(xp);
    mean_income = total_income / total_pop;

    % -------------------------------
    % Step 4: cumulative income
    % -------------------------------
    cum_xp = cumsum(xp);

    % -------------------------------
    % Step 5: Gini coefficient
    % -------------------------------
    weighted_area = 0.0;
    for i = 2:N
        weighted_area = weighted_area + (cum_xp(i-1) + cum_xp(i)) * pop(idx(i));
    end

    G = 1.0 - weighted_area / (mean_income * total_pop * total_pop);

end