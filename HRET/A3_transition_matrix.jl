# Module for distribution

module transition_matrix

export BTM, BTM_housing

function TM(a_grid, a_prime, n)    
    # transition matrix for (a, a')
    # 초기화
    transition_matrix = zeros(n, n)

    # 각 a_prime 값에 대해 transition matrix 채우기
    for i = 1:length(a_prime)
        # a_prime 값이 a_grid 범위 내에 있는지 확인
        if a_prime[i] >= a_grid[1] && a_prime[i] <= a_grid[end]
            # a_prime이 a_grid 내에 있으면 해당하는 인덱스 찾기
            idx = argmin(abs.(a_grid .- a_prime[i]))

            # transition matrix 업데이트
            if idx <= n # a_prime 값이 a_grid의 범위 내에 있으면
                # a_prime이 a_grid 값과 정확히 일치할 때
                if a_prime[i] == a_grid[idx]
                    transition_matrix[i, idx] = 1
                elseif a_prime[i] > a_grid[idx] # 가까운 쪽이 왼쪽에 있는 경우
                    ratio = (a_prime[i] - a_grid[idx]) / (a_grid[idx+1] - a_grid[idx])
                    transition_matrix[i, idx] = 1 - ratio
                    transition_matrix[i, idx+1] = ratio
                elseif a_prime[i] < a_grid[idx]  # 가까운 쪽이 오른쪽에 있는 경우
                    ratio = (a_prime[i] - a_grid[idx-1]) / (a_grid[idx] - a_grid[idx-1])
                    transition_matrix[i, idx] = ratio
                    transition_matrix[i, idx-1] = 1 - ratio
                end
            else # a_prime 값이 a_grid의 최댓값보다 큰 경우
                transition_matrix[i, i] = 1
            end
        end
    end

    return transition_matrix
end;

function BTM(a_grid, pol_a, Π, N_a, N_e)
    B = zeros(N_e*N_a, N_e*N_a)
    for j=1:N_e
        block = TM(a_grid, pol_a[j,:], N_a)    
        # block indexing
        r_start = N_a*(j-1)+1
        r_end   = N_a*j
        for k=1:N_e
            c_start = N_a*(k-1)+1
            c_end   = N_a*k
            B[r_start:r_end,c_start:c_end]=Π[j,k].*block
        end
    end
    return B
end


# aggregation을 위한 transtion matrix 생성 함수
function BTM_housing(pol_a, pol_h, Pi, N,J,K, a_grid, h_grid)
    # 1단계 a, h의 transtion block 만들기
    block_transition = zeros(N*J,N*J,K)
    for k=1:K
        for n=1:N
            for j=1:J
            
                # 값 알아보기 쉽게하기
                h_prime = pol_h[n,j,k]
                a_prime = pol_a[n,j,k]
            
                # h->h' 위치찾기
                h_idx = argmin(abs.(h_grid .- h_prime))                
                # a->a' 위치찾기
                a_idx = argmin(abs.(a_grid .- a_prime))           
                row_num = (j-1)*N + n
                col_num = (h_idx-1)*N + a_idx
                
                if a_prime == a_grid[a_idx] # a'이 a 그리드값과 동일
                    block_transition[row_num, col_num, k] = 1
                
                elseif a_prime > a_grid[a_idx] # a그리드랑 가까운 쪽이 왼쪽인 경우
                    ratio = (a_prime - a_grid[a_idx]) / (a_grid[a_idx+1] - a_grid[a_idx])
                    block_transition[row_num, col_num, k] = 1 - ratio
                    block_transition[row_num, col_num+1, k] = ratio
                
                elseif a_prime < a_grid[a_idx] # a그리드랑 가까운 쪽이 오른쪽인 경우
                    ratio = (a_grid[a_idx] - a_prime) / (a_grid[a_idx] - a_grid[a_idx-1])
                    block_transition[row_num, col_num, k] = 1- ratio
                    block_transition[row_num, col_num-1, k] = ratio
                end
            
            end
        end    
    end
    
    
    # 2단계 block matrix를 배치하기
    W = zeros(N*J*K,N*J*K);

    for k1=1:K
        for k2=1:K
            # block indexing
            r_start = N*J*(k1-1)+1
            r_end   = N*J*k1
            c_start = N*J*(k2-1)+1
            c_end   = N*J*k2
            W[r_start:r_end, c_start:c_end] = Pi[k1,k2].*block_transition[:,:,k1];
        end
    end
    
    return W
end;


end