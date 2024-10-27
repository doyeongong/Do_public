# Newton-Raphson method
# edited by Yeongwoong Do(2024)

module maximize_simple

export newton_one_maximize

function newton_one_maximize(x_ini, f, x_lower, x_upper; h=1e-6, tol=1e-8, max_iter=20)
    x = copy(x_ini)
    
    for i in 1:max_iter
        try
            fx = f(x)
        catch
            return x, -1e10 
        end
        
        # 1차 도함수(수치 미분)
        fx_h = f(x + h)
        first_derivative = (fx_h - fx) / h

        # 2차 도함수(수치 미분)
        fx_2h = f(x + 2h)
        second_derivative = (fx_2h - 2fx_h + fx) / h^2

        # 수렴 조건: 1차 도함수가 tol 이하일 때 멈춤
        if abs(first_derivative) < tol
            break
        end

        # 2차 도함수가 0이거나 양수이면, 최대값이 아닌 지점이므로 멈춤
        if second_derivative >= 0
            # println("error! No maximum found (may be a minimum or saddle point).")
            break
            return x, -1e10 
        end
        
        # Newton's method로 업데이트
        x = x - first_derivative / second_derivative

        # 범위 제약을 적용
        if x < x_lower
            x = x_lower
            break
            return x, -1e10
        elseif x > x_upper
            x = x_upper
            break
            return x, -1e10
        end
    end

    # 최적화된 함수값 계산
    fx_max = f(x)

    return x, fx_max
end

end