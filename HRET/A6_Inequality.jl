# Module for HA with housimg model help functions
# Yeongwoong Do(2024)

module Inequality

export Gini, Lorenze, logper, perentile_10

# Gini 계산기
function Gini(pop, x, tol)
    # pop의 원소가 0보다 큰 인덱스만 필터링
    # 이것을 하지 않으면 아무도 비싼 자산을 안가지고 있는데, 가지고 있는 것처럼 나와서 지니계수가 이상하게 나옴
    valid_idx = findall(p -> p > tol, pop)
    # 유효한 인덱스를 사용하여 pop과 x 벡터를 필터링
    filtered_pop = pop[valid_idx]
    filtered_x = x[valid_idx]

    # x를 기준으로 정렬된 인덱스를 얻기
    idx = sortperm(filtered_x)
    # 인덱스를 사용하여 두 벡터 모두 정렬
    sorted_x = filtered_x[idx]
    sorted_pop = filtered_pop[idx]
    
    # gini 계수 시산
    N = length(filtered_pop)
    sorted_x = sorted_x./sum(sorted_x)
    sorted_pop = sorted_pop./sum(sorted_pop)
    cum_pop = zeros(N)
    cum_x = zeros(N)
    cum_pop[1] = sorted_pop[1]
    cum_x[1] = sorted_x[1]
    under_area = sorted_x[1]*sorted_pop[1]*0.5
    for i=2:N
        cum_pop[i]=cum_pop[i-1] + sorted_pop[i]
        cum_x[i]=cum_x[i-1] + sorted_x[i]
        under_area = under_area + 0.5*sorted_pop[i]*sorted_x[i] + cum_x[i-1]*sorted_pop[i]
    end
    gini = (0.5-under_area)/0.5
    return gini
end;


function Lorenze(pop,x, tol)
    # pop의 원소가 0보다 큰 인덱스만 필터링
    # 이것을 하지 않으면 아무도 비싼 자산을 안가지고 있는데, 가지고 있는 것처럼 나와서 지니계수가 이상하게 나옴
    valid_idx = findall(p -> p > tol, pop)
    # 유효한 인덱스를 사용하여 pop과 x 벡터를 필터링
    filtered_pop = pop[valid_idx]
    filtered_x = x[valid_idx]

    # x를 기준으로 정렬된 인덱스를 얻기
    idx = sortperm(filtered_x)
    # 인덱스를 사용하여 두 벡터 모두 정렬
    sorted_x = filtered_x[idx]
    sorted_pop = filtered_pop[idx]
    
    # gini 계수 시산
    N = length(filtered_pop)
    sorted_x = sorted_x./sum(sorted_x)
    sorted_pop = sorted_pop./sum(sorted_pop)
    cum_pop = zeros(N)
    cum_x = zeros(N)
    cum_pop[1] = sorted_pop[1]
    cum_x[1] = sorted_x[1]
    under_area = sorted_x[1]*sorted_pop[1]*0.5
    for i=2:N
        cum_pop[i]=cum_pop[i-1] + sorted_pop[i]
        cum_x[i]=cum_x[i-1] + sorted_x[i]
        under_area = under_area + 0.5*sorted_pop[i]*sorted_x[i] + cum_x[i-1]*sorted_pop[i]
    end
    gini = (0.5-under_area)/0.5
    return gini, cum_pop, cum_x
end;

# p9010 계산기
function logper(pop, x, tol)
    valid_idx = findall(p -> p > tol, pop)
    filtered_pop = pop[valid_idx]
    filtered_x = x[valid_idx]

    # x를 기준으로 정렬된 인덱스를 얻기
    idx = sortperm(filtered_x)
    # 인덱스를 사용하여 두 벡터 모두 정렬
    sorted_x = filtered_x[idx]
    sorted_pop = filtered_pop[idx]

    cum_pop = cumsum(sorted_pop)./sum(sorted_pop)
    grid = [0.99 0.90 0.50 0.10]
    p = zeros(length(grid))
    
    for (idx, val) in enumerate(grid)
        # 타겟 값
        target = val
        # 각 원소와 타겟 값의 차이의 절대값을 계산
        differences = abs.(cum_pop .- target)
        # 가장 작은 차이와 그 위치 찾기
        ~, min_idx = findmin(differences)
        p[idx] = log(sorted_x[min_idx])
    end

    return p[1], p[2], p[3], p[4]
end;

# percentile 계산기
function perentile_10(pop, x, tol)
    valid_idx = findall(p -> p > tol, pop)
    filtered_pop = pop[valid_idx]
    filtered_x = x[valid_idx]

    # x를 기준으로 정렬된 인덱스를 얻기
    idx = sortperm(filtered_x)
    # 인덱스를 사용하여 두 벡터 모두 정렬
    sorted_x = filtered_x[idx]
    sorted_pop = filtered_pop[idx]

    cum_pop = cumsum(sorted_pop)./sum(sorted_pop)
    grid = range(0.10, 0.90, 9)
    p = zeros(length(grid))
    
    for (idx, val) in enumerate(grid)
        # 타겟 값
        target = val
        # 각 원소와 타겟 값의 차이의 절대값을 계산
        differences = abs.(cum_pop .- target)
        # 가장 작은 차이와 그 위치 찾기
        ~, min_idx = findmin(differences)
        p[idx] = sorted_x[min_idx]
    end

    return p
end;


end