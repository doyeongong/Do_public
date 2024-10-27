# Newton-Raphson method
# edited by Yeongwoong Do(2024)

module root_finding

export newton_one_root, newton_two_root, newton_one_root_secant, bisection_one_root

function newton_one_root(x_ini, f; h=1e-6, tol=1e-8, max_iter=20)

    x = copy(x_ini)
    
    for i in 1:max_iter
        fx = f(x)
        if abs(fx) < tol
            break
        end
        fx_h = f(x + h)
        derivative = (fx_h - fx) / h
        
        if derivative == 0
            println("Derivative is zero, cannot proceed further.")
            break
        end

        x = x - (fx) / derivative
    end

    return x
end

function newton_two_root(x_ini, f; k=0.5, h=1e-6, tol=1e-8, max_iter=20)
    
    x  = copy(x_ini)
    y  = copy(x_ini) 
    
    for i in 1:max_iter
        y[1], y[2] = f(x[1], x[2])
        diff = maximum(abs.(y))
        
        if diff < tol
            println("Converged to target within tolerance.")
            break
        end

        y_11, y_21 = f(x[1]+h, x[2])
        y_12, y_22 = f(x[1], x[2]+h)
        J = [(y_11-y[1])/h (y_12-y[1])/h; 
              (y_21-y[2])/h (y_22-y[2])/h]     

        x = x - k*inv(J)*vec(y)
    end

    println("solution $(x)")
    return x
end

function newton_one_root_secant(x0, x1, f; tol=1e-8, max_iter=20)
    # 초기값 설정
    x_prev = copy(x0)
    x = copy(x1)
    fx_prev = f(x_prev)
    fx = f(x)

    for i in 1:max_iter
        # 수렴 여부 체크
        if abs(fx) < tol
            break
        end
        
        # 평균 변화율(할선법)을 이용한 도함수 근사
        derivative = (fx - fx_prev) / (x - x_prev)
        
        if derivative == 0
            println("Derivative is zero, cannot proceed further.")
            break
        end

        # 둘 중 0에 가까운 점을 사용
        if abs(fx_prev) < abs(fx)
            x_can, fx_can = x_prev, fx_prev
        else
            x_can, fx_can = x, fx
        end

        # Newton 방법과 유사하게 x를 업데이트
        x_new = x_can - fx_can / derivative
        
        # 값 업데이트
        x_prev, fx_prev = x, fx
        x, fx = x_new, f(x_new)
    end

    return x
end

function bisection_one_root(x_ini, f; tol=1e-8, max_iter=30)

    x = copy(x_ini)
    
    for i in 1:max_iter
        x_in = (x[1]+x[2])/2
        fx = f(x_in)

        if abs(fx) < tol
            break
            return x_in
        elseif fx>0
            x[1] = x_in
        elseif fx<0
            x[2] = x_in
        end
    end

    return (x[1]+x[2])/2
end

end