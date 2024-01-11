PROGRAM Partial_retire
  
    IMPLICIT NONE
    INTEGER :: i, j, k
    
    ! Parameter
    REAL(4), PARAMETER :: sigma = 2.000d0  		    ! risk aversion
    REAL(4), PARAMETER :: beta = 0.9750d0   		    ! discount factor
    REAL(4), PARAMETER :: alpha = 0.7780d0  		    ! consumption share
    REAL(4), PARAMETER :: theta = 0.400d0  		    ! minimum down payment(1-LTV)
    REAL(4), PARAMETER :: tauselling = 0.0700d0    ! selling cost
    REAL(4), PARAMETER :: taubuying = 0.0250d0     ! buying cost
    REAL(4), PARAMETER :: r = 0.02370d0      		    ! risk-free rate
    REAL(4), PARAMETER :: kappa = 0.01330d0  		    ! spread
    REAL(4), PARAMETER :: delta0 = 0.015000d0      	! depreciation of owning house
    REAL(4), PARAMETER :: deltaR = 0.02000d0      	! depreciation of rental house
    REAL(4), PARAMETER :: tauH = 0.0010d0        	! property tax
    REAL(4), PARAMETER :: tauM = 1.000d0        	! Mortgage induction
    REAL(4), PARAMETER :: tauC = 0.0000d0         ! Comprehensive real estate tax
    REAL(4), PARAMETER :: phi = 0.30000d0           ! landlord fixed cost
    REAL(4), PARAMETER :: tauY_adj = 1.0250000d0           ! tau_Y adjustment under Gov't BC holding case
  
    ! Heathcote(2005) calibration
    INTEGER, PARAMETER :: numwage = 7
    INTEGER, PARAMETER :: numhousing = 7
    INTEGER, PARAMETER :: numequity = 60
    INTEGER, PARAMETER :: numshelter = 14
    
    REAL, DIMENSION(numwage) :: wage, tauY
    REAL, DIMENSION(numhousing) :: housing
    REAL, DIMENSION(numshelter) :: shelter
    REAL, DIMENSION(1:numequity) :: equity
    REAL, DIMENSION(numwage, numwage) :: wage_trans
    REAL, DIMENSION(numwage, numhousing, numequity) :: v_ini=0.00d0, v_old, v_next
    INTEGER, DIMENSION(numwage, numhousing, numequity) :: h_prime_index, dm_prime_index, shelter_index
    
    ! VFI 루프 필요한 변수들 선언
    REAL :: error=10.0d0, vtol=1e-6
    integer(4) :: counter=0, max_iter=1200, w, h, dm, h_prime, s, dm_prime
    integer(4) :: zeroequityind, min_dm_prime
    real :: borr, maxval_dm_prime
    Integer(4) :: maxindex_h_prime, supindex_dm_prime, supindex_shelter
    real :: maxval, value, maxval_shelter, income_tax=0.00d0
    Integer(4) :: maxindex_shelter, maxindex_dm_prime, max2index_dm_prime
    real :: consumption, Idummy, M, ytilde, payroll_taxes
    real :: deposit, mortgage_debt, income, property_tax=0.0d0, cret = 0.0d0, utility, change, expval
    real :: landlord_dummy=0.00d0, ownhouse_dummy=0.00d0
    real :: global_min=-99999999999.99
    real :: start_time, end_time, elapsed_time
  
    ! Aggregation simulation에 필요한 변수들 선언
    INTEGER, PARAMETER :: population=10000, time=5000, burn_in=2000
    INTEGER, DIMENSION(population, time) :: w_i_agent=0, h_i_agent=0, dm_i_agent=0, s_i_agent=0
    REAL, DIMENSION(population, time) :: w_agent, h_agent, dm_agent, s_agent
    REAL, DIMENSION(time) :: sum_wage=0.0d0, sum_h=0.0d0, sum_dm=0.0d0, sum_s=0.0d0, sum_h4=0.0d0, sum_s4=0.0d0
    REAL, DIMENSION(time) :: sum_debt=0.0d0, sum_deposit=0.0d0, sum_owner=0.0d0, sum_landlord=0.0d0
    REAL, DIMENSION(time) :: sum_trade= 0.0d0, sum_mi = 0.0d0
    REAL :: random_numbers(population)
    REAL, DIMENSION(numwage,numwage) :: cum_trans
    INTEGER :: w_colum, w_row, t, seed_num
    REAL :: avg_W=0, H_d=0, S_d=0, Debt=0, Saving=0
    REAL :: share_homeowners=0, share_renters=0, share_landlords=0
    REAL :: share_h4=0, share_s4=0, share_trading=0, share_mortgage=0
  
    ! bisection-grid
    REAL(4) :: rho_lower=0.145000000d0, rho, q_lower=1.3610d0, q
    INTEGER(4), PARAMETER :: rho_grid=1, q_grid=1
    REAL, DIMENSION(rho_grid) :: rho_vec
    REAL, DIMENSION(q_grid) :: q_vec
    CHARACTER(LEN=30) :: filename
    
    filename = 'grid_7_new.txt'
    ! 파일 열기
    OPEN(UNIT=10, FILE=filename, STATUS='UNKNOWN')
  
    ! start time
    call CPU_TIME(start_time)
  
    ! start time
    call CPU_TIME(start_time)
    rho_vec(1) = rho_lower
    q_vec(1) = q_lower
    do j=2,rho_grid
      rho_vec(j) = rho_vec(j-1) + 0.0010d0
    enddo
    do i=2,q_grid
      q_vec(i) = q_vec(i-1) + 0.0050d0
    enddo
  
    print *, "------------start calculation------------"
     ! input
    wage(1) = 0.12720d0  
    wage(2) = 0.37450d0
    wage(3) = 0.67730d0  
    wage(4) = 1.00000d0
    wage(5) = 1.55550d0
    wage(6) = 1.85190d0  
    wage(7) = 2.78880d0  
      
    housing(1) = 0.00d0
    housing(2) = 1.00d0  
    housing(3) = 2.00d0
    housing(4) = 3.00d0
    housing(5) = 5.00d0
    housing(6) = 7.00d0
    housing(7) = 10.00d0
    
    shelter(1) = 0.25d0
    shelter(2) = 0.50d0
    shelter(3) = 0.750d0
    shelter(4) = 1.00d0
    shelter(5) = 1.250d0
    shelter(6) = 1.500d0
    shelter(7) = 1.750d0
    shelter(8) = 2.00d0
    shelter(9) = 2.500d0
    shelter(10) = 3.00d0
    shelter(11) = 4.00d0
    shelter(12) = 5.00d0
    shelter(13) = 7.00d0
    shelter(14) = 10.00d0
  
    wage_trans(1,1) = 0.4140d0
    wage_trans(1,2) = 0.5751d0
    wage_trans(1,3) = 0.0109d0
    wage_trans(1,4) = 0.0000d0
    wage_trans(1,5) = 0.0000d0
    wage_trans(1,6) = 0.0000d0
    wage_trans(1,7) = 0.0000d0
    
    wage_trans(2,1) = 0.1056d0
    wage_trans(2,2) = 0.4646d0
    wage_trans(2,3) = 0.4199d0
    wage_trans(2,4) = 0.0098d0
    wage_trans(2,5) = 0.0000d0
    wage_trans(2,6) = 0.0000d0
    wage_trans(2,7) = 0.0000d0
    
    wage_trans(3,1) = 0.0034d0
    wage_trans(3,2) = 0.1395d0
    wage_trans(3,3) = 0.5187d0
    wage_trans(3,4) = 0.3308d0
    wage_trans(3,5) = 0.0075d0
    wage_trans(3,6) = 0.0000d0
    wage_trans(3,7) = 0.0000d0
    
    wage_trans(4,1) = 0.0000d0
    wage_trans(4,2) = 0.0052d0
    wage_trans(4,3) = 0.1732d0
    wage_trans(4,4) = 0.5339d0
    wage_trans(4,5) = 0.2825d0
    wage_trans(4,6) = 0.0052d0
    wage_trans(4,7) = 0.0000d0
    
    wage_trans(5,1) = 0.0000d0
    wage_trans(5,2) = 0.0000d0
    wage_trans(5,3) = 0.0075d0
    wage_trans(5,4) = 0.2028d0
    wage_trans(5,5) = 0.5187d0
    wage_trans(5,6) = 0.2676d0
    wage_trans(5,7) = 0.0034d0
  
    wage_trans(6,1) = 0.0000d0
    wage_trans(6,2) = 0.0000d0
    wage_trans(6,3) = 0.0000d0
    wage_trans(6,4) = 0.0098d0
    wage_trans(6,5) = 0.2186d0
    wage_trans(6,6) = 0.4646d0
    wage_trans(6,7) = 0.3068d0
  
    wage_trans(7,1) = 0.3778d0
    wage_trans(7,2) = 0.0000d0
    wage_trans(7,3) = 0.0000d0
    wage_trans(7,4) = 0.0000d0
    wage_trans(7,5) = 0.0109d0
    wage_trans(7,6) = 0.1974d0
    wage_trans(7,7) = 0.4140d0
  
    tauY(1) = 0.020d0 * tauY_adj
    tauY(2) = 0.032d0 * tauY_adj
    tauY(3) = 0.058d0 * tauY_adj
    tauY(4) = 0.080d0 * tauY_adj
    tauY(5) = 0.1080d0 * tauY_adj
    tauY(6) = 0.1270d0 * tauY_adj
    tauY(7) = 0.1780d0 * tauY_adj
        
    equity(1) = -15.00d0
    
    DO i = 2, 29
      equity(i) = equity(i - 1) + 0.50d0
    END DO 
    DO i = 20, 49
      equity(i) = equity(i - 1) + 0.1d0
    END DO
    DO i = 50, 60
      equity(i) = equity(i - 1) + 0.50d0
    END DO
    equity(39)=0.00d0
    zeroequityind=39
  
    Do j = 1, q_grid
      q = q_vec(j)
    DO k = 1, rho_grid
      rho = rho_vec(k)
      print *, "=====New iteration======"
      print *, "houseprice =", q, "rent = ", rho
      
      ! zero로 초기화
      avg_W=0.0d0
      H_d=0.0d0
      S_d=0.0d00
      Debt=0.0d0
      Saving=0.0d0
      sum_wage=0.0d0
      sum_h=0.0d0
      sum_s=0.0d0
      sum_debt=0.0d0
      sum_deposit=0.0d0
      sum_owner=0.0d0
      sum_landlord=0.0d0
      sum_h4=0.0d0
      sum_s4=0.0d0
      sum_trade = 0.0d0
      sum_mi  = 0.0d0
  
    counter=0
    error = 10.0
    v_old = v_ini
  
  
    v_old = v_ini
  
    DO while (counter < max_iter .AND. error > vtol)
      change=0.00d0
      Do w=1,numwage
      Do h=1,numhousing
      Do dm=1,numequity
          maxval = global_min
          maxindex_h_prime = 1
          supindex_dm_prime = zeroequityind
          supindex_shelter = 1
          
          Do h_prime=1,numhousing
        maxval_shelter = global_min
        maxindex_shelter = 1
        max2index_dm_prime = zeroequityind
        
        do s = 1, numshelter
          maxval_dm_prime = global_min
          maxindex_dm_prime = zeroequityind
  
          do dm_prime = 1, numequity
  
            if (housing(h_prime) .gt. 0.0 .and. shelter(s) .gt. housing(h_prime)) then
              value = global_min
            elseif (equity(dm_prime) .lt. -(1-theta)*q*housing(h_prime)) then
              value = global_min
            else
              ! 주택거래 발생여부 판별
              if (h_prime .ne. h) then ! 주택거래발생 더미
                Idummy = 1.00d0
              else ! 거래 없음
                Idummy = 0.00d0
              endif
  
              ! 은행 포지션 판별
              if (equity(dm) .le. 0.00d0) then
                ! vs 은행 포지션이 차입자인 경우
                deposit = 0.0d0
                mortgage_debt = -equity(dm) 
              else
                ! vs 은행 포지션이 대부자인 경우
                deposit = equity(dm)
                mortgage_debt = 0.0d0
              endif
  
              ! 임대업자
              if (housing(h_prime) .gt. shelter(s)) then
                landlord_dummy=1.00d0
              else
                landlord_dummy=0.00d0
              endif
  
              ! 주택보유
              if (housing(h_prime) .gt. 0.0) then
                ownhouse_dummy=1.00d0
              else
                ownhouse_dummy = 0.00d0
              endif
  
              ! tax 처리
              income = wage(w) + r*deposit + landlord_dummy*rho*(housing(h_prime)-shelter(s))
              property_tax = tauH*q*housing(h_prime)
              if (h_prime .eq. numhousing) then
                  cret = tauC*q*housing(h_prime)
              else
                  cret = 0.00d0
              endif
              call tax_induction(w, income, housing(h_prime), shelter(s), mortgage_debt, r+kappa, &
              landlord_dummy, property_tax, tauM, deltaR, q, income_tax, tauY)
  
              ! 유지비용
              M = delta0*shelter(s) + deltaR*MAX(housing(h_prime)-shelter(s), 0.0d0)
  
              ! 소비계산
              consumption = wage(w) + (1.00d0 + r)*deposit + q*housing(h) &
              -(1.00d0+ r +kappa)*mortgage_debt &
              - rho*shelter(s) - (q-rho)*housing(h_prime) - equity(dm_prime) &
              - Idummy*tauselling*q*housing(h) &
              - Idummy*taubuying*q*housing(h_prime) &
              - property_tax -cret - income_tax - ownhouse_dummy*q*M - landlord_dummy*phi
              
              if (consumption .lt. 0.0) then
                value = global_min
              else
                call utilityflow(consumption, shelter(s), utility, alpha, sigma)
                call expected_val(w, h_prime, dm_prime, v_old, wage_trans, numwage, numhousing, numequity, zeroequityind, expval)
                value = utility + beta* expval
              end if
  
            end if
  
            if (value .gt. maxval_dm_prime) then
              maxval_dm_prime = value
              maxindex_dm_prime = dm_prime
            endif
          enddo
        
          if (maxval_dm_prime .gt. maxval_shelter) then
            maxval_shelter = maxval_dm_prime
            maxindex_shelter = s
            max2index_dm_prime = maxindex_dm_prime
          endif    
            enddo
  
        if (maxval_shelter .gt. maxval) then
          maxval = maxval_shelter
          maxindex_h_prime = h_prime
          supindex_shelter = maxindex_shelter
          supindex_dm_prime = max2index_dm_prime
        endif
  
        enddo
  
      change = change + ABS(maxval - v_old(w,h,dm))
  
      v_next(w,h,dm) = maxval
      h_prime_index(w,h,dm) = maxindex_h_prime
      dm_prime_index(w,h,dm) = supindex_dm_prime
      shelter_index(w,h,dm) = supindex_shelter
  
    enddo
    enddo
    enddo
    error = change
    counter= counter +1
    v_old = v_next
  enddo
  
  write(*,*) "iteration :", counter
  write(*,*) "error value: ", error
  write(*,*) "The richest less saving?", dm_prime_index(numwage,numhousing,numequity) .lt. numequity
  print *, "----------------Aggregation---------------------"
  
  ! Aggregation simulation
  ! initial value for a simulation
  do i=1,population
      w_i_agent(i,1)=2
      h_i_agent(i,1)=1
      dm_i_agent(i,1)=zeroequityind
      s_i_agent(i,1)=shelter_index(w_i_agent(i,1),h_i_agent(i,1), dm_i_agent(i,1))
      w_agent(i,1)=wage(w_i_agent(i,1))
      h_agent(i,1)=housing(h_i_agent(i,1))
      dm_agent(i,1)=equity(dm_i_agent(i,1))
      s_agent(i,1) = shelter(s_i_agent(i,1))
      ! aggregation
      sum_wage(1) = sum_wage(1) + w_agent(i,1)
      sum_h(1) = sum_h(1) + h_agent(i,1)
      sum_s(1) = sum_s(1) + s_agent(i,1)
  enddo
  
  ! cumulative wage transition matrix
  do w_row = 1, numwage
      do w_colum = 1, numwage
          if (w_colum .eq. 1) then
              cum_trans(w_row,w_colum)=wage_trans(w_row,w_colum)
          else
              cum_trans(w_row,w_colum)=cum_trans(w_row,w_colum-1) + wage_trans(w_row,w_colum)
          endif
  enddo
  enddo
  
  ! simulation
  seed_num = 20230725
  call random_seed(seed_num)
  
  print *, "----initialized completed----"
  
  do t=2,time
  
      ! 난수 생성
      call RANDOM_NUMBER(random_numbers)
  
      do i=1,population
          ! wage 
          call prob_match(numwage,random_numbers(i), cum_trans(w_i_agent(i,t-1),:),w_i_agent(i,t))
          if ((w_i_agent(i,t-1) .eq. numwage) .AND. (w_i_agent(i,t) .eq. 1)) then ! 새로 태어난 경우
            h_i_agent(i,t)=1  ! 집은 없고
            dm_i_agent(i,t)=zeroequityind ! 돈도 없다
          else ! 나머지는 그대로
            h_i_agent(i,t)=h_prime_index(w_i_agent(i,t-1), h_i_agent(i,t-1),dm_i_agent(i,t-1))
            dm_i_agent(i,t)=dm_prime_index(w_i_agent(i,t-1), h_i_agent(i,t-1),dm_i_agent(i,t-1))
          end if
          s_i_agent(i,t)=shelter_index(w_i_agent(i,t-1), h_i_agent(i,t-1),dm_i_agent(i,t-1))   
          ! value
          w_agent(i,t)=wage(w_i_agent(i,t))
          h_agent(i,t)=housing(h_i_agent(i,t))
          dm_agent(i,t)=equity(dm_i_agent(i,t))
          s_agent(i,t) = shelter(s_i_agent(i,t))
          ! Aggregation
          sum_wage(t) = sum_wage(t) + w_agent(i,t)
          sum_h(t) = sum_h(t) + h_agent(i,t)
          sum_s(t) = sum_s(t) + s_agent(i,t)
          if (h_agent(i,t) .ne. h_agent(i,t-1)) then
            if ((w_i_agent(i,t-1) .eq. numwage) .AND. (w_i_agent(i,t) .eq. 1)) then
              sum_trade(t) = sum_trade(t)
            else
              sum_trade(t) = sum_trade(t) + 1
            endif
          endif
          if (dm_agent(i,t) .lt. 0.0d0) then
              sum_debt(t) = sum_debt(t) + dm_agent(i,t)
              sum_mi(t) = sum_mi(t) + 1
          else
              sum_deposit(t) = sum_deposit(t) + dm_agent(i,t)
          endif
          if (h_agent(i,t) .gt. 0.0) then
            sum_owner(t) = sum_owner(t) + 1
          end if
          if (h_agent(i,t) .gt. s_agent(i,t)) then
            sum_landlord(t) = sum_landlord(t) + 1
          end if
          if (h_i_agent(i,t) .eq. numhousing) then
            sum_h4(t) = sum_h4(t) + 1
          end if
          if (s_i_agent(i,t) .eq. numshelter) then
            sum_s4(t) = sum_s4(t) + 1
          end if
      enddo
  enddo
  print *, "----end of simulation----"
  
  avg_W = SUM(sum_wage(burn_in+1:time))/(population*(time-burn_in))
  H_d = SUM(sum_h(burn_in+1:time))/(population*(time-burn_in))
  S_d = SUM(sum_s(burn_in+1:time))/(population*(time-burn_in))
  Debt = SUM(sum_debt(burn_in+1:time))/(population*(time-burn_in))
  Saving = SUM(sum_deposit(burn_in+1:time))/(population*(time-burn_in))
  share_homeowners = SUM(sum_owner(burn_in+1:time))/(population*(time-burn_in))
  share_renters = 1 - share_homeowners
  share_h4 = SUM(sum_h4(burn_in+1:time))/(population*(time-burn_in))
  share_s4 = SUM(sum_s4(burn_in+1:time))/(population*(time-burn_in))
  share_landlords = SUM(sum_landlord(burn_in+1:time))/(population*(time-burn_in))
  share_trading = SUM(sum_trade(burn_in+1:time))/(population*(time-burn_in))
  share_mortgage = SUM(sum_mi(burn_in+1:time))/(population*(time-burn_in))
  
  write(*,*) "Average wage :", avg_w
  write(*,*) "Aggregate housing demand :", H_d
  write(*,*) "Aggregate Shelter demand :", S_d
  write(*,*) "Aggregate debt :", Debt
  write(*,*) "Aggregate deposit :", Saving
  write(*,*) "share of renters :", share_renters
  write(*,*) "share of landlords :", share_landlords
  write(*,*) "share of CRET payers :", share_h4
  write(*,*) "share of high value shelter :", share_s4
  write(*,*) "share of trading :", share_trading
  write(*,*) "share of mortgage indebted :", share_mortgage
  
  WRITE(10, '(F6.4)') q, rho, H_d, S_d, share_Renters, share_landlords, share_h4, share_trading, share_mortgage
  enddo
  enddo
  
  CLOSE(10)
  
  ! ending time
  call CPU_TIME(end_time)
  elapsed_time = end_time - start_time
  print *, "Elapsed time: ", elapsed_time, "seconds"
  
  
  END PROGRAM Partial_retire
  
  
    subroutine tax_induction(w, y, h_prime, s, mortgage, mortgage_rate, landlord_dummy, &
       property_tax, tauM, deltaR, q, income_tax, tauY)
      ! note: 여기서 h_prime과 s는 index가 아닌 값이 들어가야 한다!
      ! 결과값 ytilde는 과세표준을 의미
      implicit none
      integer :: w
      real :: y, h_prime, s, mortgage, mortgage_rate, mortgage_dummy, landlord_dummy
      real :: property_tax, tauM, deltaR, q, rho, stand_ded, ytilde, income_tax, deduction
      real, DIMENSION(7) :: tauY
    
      stand_ded = 0.0690d0    ! 기본공제 300만원
     
      mortgage_dummy = 0.00d0
  
      if (h_prime .le. 7.00d0 ) then
        mortgage_dummy = 1.00d0
      end if
    
      ytilde = y + mortgage_dummy*(-tauM*mortgage_rate*mortgage) - stand_ded
    
      if (ytilde .lt. 0.0) then ! 음수라면 0을 return
        ytilde = 0.0d0
      endif
    
      if (w .eq. 1) then
        income_tax = max(tauY(1)*ytilde,0d0)
      elseif (w .eq. 2) then
        income_tax = max(tauY(2)*ytilde,0d0)
      elseif (w .eq. 3) then
        income_tax = max(tauY(3)*ytilde,0d0)
      elseif (w .eq. 4) then
        income_tax = max(tauY(4)*ytilde,0d0)
      elseif (w .eq. 5) then 
        income_tax = max(tauY(5)*ytilde,0d0)
      elseif (w .eq. 6) then 
        income_tax = max(tauY(6)*ytilde,0d0)
      elseif (w .eq. 7) then 
        income_tax = max(tauY(7)*ytilde,0d0)
      endif
    
      return
      end
    
    subroutine utilityflow(c, s, utility, alpha, sigma)
      implicit none
      real :: c, s, utility, alpha, sigma
      utility = (c**(alpha)*s**(1-alpha))**(1-sigma)/(1-sigma)
      return
      end
    
    subroutine expected_val(w, h_prime, dm_prime, v_guess, wage_trans, &
      numwage, numhousing, numequity, zeroequityind, expval)
      implicit none
      integer :: w, w_prime, h_prime, dm_prime, numwage, numhousing, numequity, zeroequityind, h_next, dm_next
      real :: expval
      REAL, DIMENSION(numwage, numhousing, numequity) :: v_guess
      REAL, DIMENSION(numwage, numwage) :: wage_trans
    
      expval = 0.00d0
    
      do w_prime=1, numwage
        if (w .le. (numwage-1)) then ! 다음기에 안죽는 경우
          expval = expval + wage_trans(w,w_prime)*v_guess(w_prime, h_prime, dm_prime)
        else 
          if (w_prime .eq. 1) then ! 다음기에 죽는 경우
            h_next = 1
            dm_next = zeroequityind
          else
            h_next = h_prime
            dm_next = dm_prime
          endif
          expval = expval + wage_trans(w,w_prime)*v_guess(w_prime, h_next, dm_next)
        endif
      enddo
    
      return
      end
    
    subroutine prob_match(dimention, var, var_grid, gridpoint)
      ! simulation 할 때, 누적분포에 맞는 값을 주는 함수
      implicit none
      integer i, dimention, gridpoint
      real var, var_grid(dimention), lower, upper
    
      if (dimention .eq. 1) then
        gridpoint = 1
      endif
    
      do i=1, dimention
        if (i .eq. 1 .AND. var .le. var_grid(i)) then
          gridpoint = 1
        elseif (i .ge. 2 .AND. (var .ge. var_grid(i-1)) .AND. (var .le. var_grid(i))) then
          gridpoint = i
        endif
        if ( var .gt. var_grid(dimention) ) then
          gridpoint = dimention
        endif
        if ( var .lt. var_grid(1) ) then
          gridpoint = 1
        endif
      enddo
      return
      end
  
  
    
  
  