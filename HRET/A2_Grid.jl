# Grid

module Grid

export log_type_grid, log_type_grid_strong_curve, densezero

function log_type_grid(g_min, g_max, n)
    pivot = abs(g_min) + 0.25 # 0에 로그를 못취하므로 -inf 방지
    grid = exp.(LinRange(log(g_min+pivot), log(g_max+pivot), n)).-pivot
    grid[1] = g_min
    return grid
end

function log_type_grid_strong_curve(g_min, g_max, n; curve_factor=2)
    pivot = abs(g_min) + 0.25 # 0에 로그를 못취하므로 -inf 방지
    grid = exp.(LinRange(log(g_min+pivot) * curve_factor, log(g_max+pivot), n)).-pivot
    grid[1] = g_min
    return grid
end


function densezero(amin, amax, N)  
    # a_grid 생성
    a_debt_grid = -(LinRange(sqrt(-amin), sqrt(0), Int(N/2))).^2
    a_delta = a_debt_grid[end] - a_debt_grid[end-1]
    a_saving_grid = (LinRange(sqrt(a_delta), sqrt(amax), Int(N/2))).^2
    a_grid = [a_debt_grid; a_saving_grid]   
    return a_grid
end

    
end