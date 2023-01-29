clear 

cd "C:\Users\82106\Documents\2. 2022학년도\가계재무건전성 평가\가계금융복지조사\가구마스터_20221012_23023_데이터"
* 이 do 파일은 T를 calibration 하는 부분을 실행하는 파일입니다.
* 이 파일은 basic_statistic.do파일로 basic_statistics.dta가 생성되어 있어야 실행이 가능합니다.
* Note: 2021년 기준 가계금융복지조사 data는 소득, 지출이 모두 2020년 중 기준
*       따라서 calibration 타겟은 2020년 은행의 가계대출 연체율 수치인 0.234%이다.

use basic_statistics.dta

* 패널 셋팅
xtset md제공용_가구고유번호 year

* 부채가구 한정
keep if debt==1
* 유동성 (2021년 자료에는 2020년 유동성이 필요)
gen LA = L.자산_금융자산
* 2020년 기준 data 사용
keep if year==2021
* 은행대출
gen 은행담보대출=부채_금융부채_담보대출_대출기관_은행금액2010년은농수협중
replace 은행담보대출=부채_금융부채_담보대출_대출기관_은행금액 if 은행담보대출==.
gen 은행대출=은행담보대출+부채_금융부채_신용대출_대출기관_은행금액
gen 은행여부=(은행대출>0)
gen 부채가중=가중값*부채
gen 은행대출가중=가중값*은행대출


******* <표 6> T값 캘리브레이션 결과 ************

* PD시산()
gen PD = .

* T=0
gen PD0 = 0 if FM>=0
replace PD0 = 0 if 0 <= LA & FM<0
replace PD0 = 1 if 0==LA & FM<0
* T=0.5
gen PD0_5 =0 if FM>=0
replace PD0_5 = 0 if (abs(FM)/12)*0.5<=LA & FM<0
replace PD0_5 = 1 - LA/((abs(FM)/12)*0.5) if (abs(FM)/12)*0.5>LA & FM<0
* T=0.69
gen PD_star =0 if FM>=0
replace PD_star = 0 if (abs(FM)/12)*0.69<=LA & FM<0
replace PD_star = 1 - LA/((abs(FM)/12)*0.69) if (abs(FM)/12)*0.69>LA & FM<0
* T=1, 3, 6, 12
foreach i in 1 3 6 12 {
	gen PD`i' =0 if FM>=0
	replace PD`i' = 0 if (abs(FM)/12)*`i'<=LA & FM<0
	replace PD`i' = 1 - LA/((abs(FM)/12)*`i') if (abs(FM)/12)*`i'>LA & FM<0
}

* 표생성
summ PD* if debt==1 [aweight=부채가중]
summ PD* if 은행여부==1 [aweight=은행대출가중]

********** <표 7> 2020년 기준 가계부문 리스크 현황  **************************

* T=0.69로 PD재시산
gen PD =0 if FM>=0
replace PD = 0 if (abs(FM)/12)*0.69<=LA & FM<0
replace PD = 1 - LA/((abs(FM)/12)*0.69) if (abs(FM)/12)*0.69>LA & FM<0
*LGD
gen LGD_ratio= PD*max((부채-자산_실물자산_부동산금액*0.8)/부채,0)
*가중치 재설정
gen 비은행대출=부채-은행대출
gen 비은행대출가중=가중값*비은행대출
gen 비은행여부=(비은행대출>0)
*소득분위
gen IQ1=(소득5분위코드보완=="Q1")
gen IQ2=(소득5분위코드보완=="Q2")
gen IQ3=(소득5분위코드보완=="Q3")
gen IQ4=(소득5분위코드보완=="Q4")
gen IQ5=(소득5분위코드보완=="Q5")

* 전체 PD EAD LGD ratio
summ PD [aweight=가중값]
summ PD LGD_ratio [aweight=부채가중]
* 은행 PD EAD LGD ratio
summ PD [aweight=가중값] if 은행여부==1
summ PD LGD_ratio [aweight=은행대출가중] if 은행여부==1
* 비은행 PD EAD LGD ratio
summ PD [aweight=가중값] if 비은행여부==1
summ PD LGD_ratio [aweight=비은행대출가중] if 비은행여부==1
* 소득분위별
forvalues i=1(1)5 {
	summ PD [aweight=가중값] if IQ`i'==1
	summ PD LGD_ratio [aweight=부채가중] if IQ`i'==1
}
