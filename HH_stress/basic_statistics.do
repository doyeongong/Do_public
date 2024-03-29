clear 
cd "C:\Users\82106\Documents\2. 2022학년도\가계재무건전성 평가\가계금융복지조사\가구마스터_20221012_23023_데이터"

***************** <표 4> 가계금융복지조사의 표본 정보****************************************
clear
import delimited "2021.csv", encoding(EUC-KR) 
save 2021_basic.dta, replace

clear
import delimited "2020.csv", encoding(EUC-KR) 
save 2020_basic.dta, replace

clear
import delimited "2019.csv", encoding(EUC-KR)
save 2019_basic.dta, replace

append using 2020_basic.dta
append using 2021_basic.dta

* 부채가구, 연체가구
gen debt=(부채>0)
gen default=(일이상원리금연체여부==1)

* 적자가구 만들기
* NI
replace 보완_처분가능소득=처분가능소득보완경상소득보완비소비지출보완 if (조사연도>=2017) & (조사연도<=2020)
drop 처분가능소득보완경상소득보완비소비지출보완
gen NI = 보완_처분가능소득
* DP
gen DP=원리금상환금액
* BLC (2021년 -> 2020년 생계비)
rename 조사연도 year
gen BLC = 632.5896 if 가구원수==1 & year==2021
replace BLC = 1077.1128 if 가구원수==2 & year==2021
replace BLC = 1393.4076 if 가구원수==3 & year==2021
replace BLC = 1709.7024 if 가구원수==4 & year==2021
replace BLC = 2021.9972 if 가구원수==5  & year==2021
replace BLC = 2342.2920 if 가구원수==6  & year==2021
replace BLC = 2660.2980 if 가구원수==7  & year==2021
replace BLC = 2660.2980+318.0060*(가구원수-7) if 가구원수>7
replace BLC = 614.5224 if 가구원수==1 & year==2020
replace BLC = 1046.3496 if 가구원수==2 & year==2020
replace BLC = 1353.6120 if 가구원수==3 & year==2020
replace BLC = 1660.8732 if 가구원수==4 & year==2020
replace BLC = 1968.1344 if 가구원수==5  & year==2020
replace BLC = 2275.3956 if 가구원수==6  & year==2020
replace BLC = 2582.6568 if 가구원수==7  & year==2020
replace BLC = 2582.6568+307.2612*(가구원수-7) if 가구원수>7
replace BLC = 601.9584 if 가구원수==1 & year==2019
replace BLC = 1024.9548 if 가구원수==2 & year==2019
replace BLC = 1325.9340 if 가구원수==3 & year==2019
replace BLC = 1626.9132 if 가구원수==4 & year==2019
replace BLC = 1927.8912 if 가구원수==5  & year==2019
replace BLC = 2228.8704 if 가구원수==6  & year==2019
replace BLC = 2529.8496 if 가구원수==7  & year==2019
replace BLC = 2529.8496+300.9792*(가구원수-7) if 가구원수>7
* FM
gen FM= NI-DP-BLC
* deficit표시
gen deficit=(FM<0)

** 표값 도출
* 전체가구수
table year, c(n md제공용_가구고유번호)
* 부채가구수
table year debt, c(n md제공용_가구고유번호)
* 적자가구수
table year deficit, c(n md제공용_가구고유번호)
* 연체가구수
table year default, c(n md제공용_가구고유번호)
* 부채가구 비중
table year [aweight=가중값], c(mean debt)
* 적자 및 연체가구 비중
table year if debt==1 [aweight=가중값], c(mean deficit mean default)

save basic_statistics.dta, replace

***********************<표 5> 연체유무별 기초통계량 비교 *****************************
clear
use basic_statistics.dta
keep if year==2021 & debt==1

* 금융자산(LA)
gen LA = 자산_금융자산


** 재무현황 부분
table year, c(med NI med DP med BLC med FM med LA)
table default, c(med NI med DP med BLC med FM med LA)
table year, c(med 자산 med 부채)
table default, c(med 자산 med 부채)

** 재무특성 및 인구사회적 특성 부분
* 기관별대출여부
gen 은행담보대출=부채_금융부채_담보대출_대출기관_은행금액2010년은농수협중
replace 은행담보대출=부채_금융부채_담보대출_대출기관_은행금액 if 은행담보대출==.
gen 은행대출=은행담보대출+부채_금융부채_신용대출_대출기관_은행금액
gen 은행여부=(은행대출>0)
gen 비은행여부=(부채-은행대출>0)
* 자영업자
gen 자영업자=0 if 가구주_종사상지위코드<3 
replace 자영업자=1  if 가구주_종사상지위코드>=3 & 가구주_종사상지위코드>=6
* 연령
gen 청년=(가구주_만연령<30)
gen 중년=(가구주_만연령>=30 & 가구주_만연령<60)
gen 장년=(가구주_만연령>=60)
gen 고학력=(가구주_교육정도_통합코드=="G4")
* 1인가구
gen 솔로가구=(가구원수==1)
keep year md제공용_가구고유번호 가중값 deficit 은행여부 비은행여부 자영업자 가구주_성별코드 청년 중년 장년 고학력 가구원수 솔로가구 default
summ [aweight=가중값]
summ [aweight=가중값] if default!=1
summ [aweight=가중값] if default==1

* Note: 변동금리 비중, 평균차입금리, 고금리비중은 계좌단위의 자료가 필요하여 RAS 시스템 상에서만 확인이 가능
