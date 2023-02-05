clear
cd "H:\Microdata\221019_최초자료\"

*** 부채자료와 가금복 메인자료와 matching test
*import delimited bc2021.csv
*rename (c1 c2 c8 c10 c13 c17) (year id 잔액 이자지급금액 최초대출금액 금리조건)
*collapse (sum) 이자지급급액 (sum) 잔액, by(id)
*save collapse_test.dta

*import delimited ggm2021.csv
*rename (c1 c2 c66 c127) (year id 부채 원리금상환_이자지급금액)
*merge m:1 year id using collapse_test.dta

***** 과거자료 
*import delimited bc2016.csv, clear
*rename (c1 c2 c9 c11 c14 c18) (year id 잔액 이자지급금액 최초대출금액 금리조건)


***** 가구별 변동금리 익스포져

import delimited bc2021.csv, clear
rename (c1 c2 c8 c10 c13 c17) (year id 잔액 이자지급금액 최초대출금액 금리조건)
gen 변동여부=(금리조건==1)
replace 변동여부=1 if 금리조건==.
replace 변동여부=0 if 금리조건==2
gen 변동이자=이자지급금액*변동여부
collapse (sum) 이자지급금액 (sum) 최초대출금액 (sum) 변동이자, by(id)
gen float_exposure=변동이자/이자지급금액
replace float_exposure=0 if float_exposure==.
save floating_exposure_2021.dta, replace

****** merging

import delimited ggm2021.csv, clear
rename (c1 c2 c66 c127) (year id 부채 원리금상환_이자지급금액)
merge m:1 id using floating_exposure.dta
save 2021_sim.dta, replace

drop if 부채==0

*** 가중치
rename (c3 c76 c100) (가중값 은행_담보 은행_신용)
gen 은행대출=은행_담보+은행_신용
gen 비은행대출=부채-은행대출
gen 부채가중=가중값*부채
gen 은행가중=가중값*은행대출
gen 비은행가중 = 가중값*비은행대출

** 2022년 금리충격 150bp
rename (c7 c152 c125 c130 c131 c132 c151 c166) (가구원수 NI DP 근로소득 사업소득 재산소득 이전소득 비소비지출)
gen 이자부담증가분=최초대출금액*float_exposure*0.015
gen DP_sim=DP+이자부담증가분

** 2022년 물가충격 5.2%
gen BLC_sim=699.105 if 가구원수==1
replace BLC_sim=1177.602 if 가구원수==2
replace BLC_sim=1519.735 if 가구원수==3
replace BLC_sim=1851.720 if 가구원수==4
replace BLC_sim=2147.618 if 가구원수==5
replace BLC_sim=2515.266 if 가구원수==6
replace BLC_sim=2855.603 if 가구원수==7
replace BLC_sim=2855.603+(324.742*(가구원수-7)) if 가구원수>7

** 2022년 주가충격 2021.3월말잔*(1-0.213)
rename (c38 c46 c57) (LA 주식채권펀드 부동산)
gen 주식채권펀드_sim=주식채권펀드*(1-0.213)
gen LA_sim=(LA-주식채권펀드)+주식채권펀드_sim
gen CA_sim=부동산*(1.044)

** 2022년 소득충격 (피보 5.6%)
* simulation 파일의 control시트 반드시 확인!!
gen total_PD=.
gen total_EAD=.
gen total_LGD=.
gen bank_PD=.
gen bank_EAD=.
gen bank_LGD=.
gen non_PD=.
gen non_EAD=.
gen non_LGD=.

gen 근로_sim=.
gen 사업_sim=.
gen prob=.
gen FM_sim=.
gen PD_sim=.
gen LGD_ratio_sim=.

forvalues i=1(1)1000 {
	set seed `i*7'
	replace prob=runiform()
	
	replace 근로_sim=근로소득*(1-0.613+0.019) if prob<0.0124
	replace 근로_sim=근로소득*(1-0.573+0.019) if prob>=0.0124 & prob<0.0306
	replace 근로_sim=근로소득*(1-0.528+0.019) if prob>=0.0306 & prob<0.0470
	replace 근로_sim=근로소득*(1-0.478+0.019) if prob>=0.0470 & prob<0.0681
	replace 근로_sim=근로소득*(1-0.423+0.019) if prob>=0.0681 & prob<0.0935
	replace 근로_sim=근로소득*(1-0.362+0.019) if prob>=0.0935 & prob<0.1230
	replace 근로_sim=근로소득*(1-0.295+0.019) if prob>=0.1230 & prob<0.1588
	replace 근로_sim=근로소득*(1-0.221+0.019) if prob>=0.1588 & prob<0.2097
	replace 근로_sim=근로소득*(1-0.139+0.019) if prob>=0.2097 & prob<0.2768
	replace 근로_sim=근로소득*(1-0.049+0.019) if prob>=0.2768 & prob<0.3765
	replace 근로_sim=근로소득*(1+0.051+0.019) if prob>=0.3765 & prob<0.5584
	replace 근로_sim=근로소득*(1+0.162+0.019) if prob>=0.5584 & prob<0.7067
	replace 근로_sim=근로소득*(1+0.284+0.019) if prob>=0.7067 & prob<0.7916
	replace 근로_sim=근로소득*(1+0.419+0.019) if prob>=0.7916 & prob<0.8487
	replace 근로_sim=근로소득*(1+0.568+0.019) if prob>=0.8487 & prob<0.8899
	replace 근로_sim=근로소득*(1+0.733+0.019) if prob>=0.8899 & prob<0.9228
	replace 근로_sim=근로소득*(1+0.916+0.019) if prob>=0.9228 & prob<0.9503
	replace 근로_sim=근로소득*(1+1.117+0.019) if prob>=0.9503 & prob<0.9746
	replace 근로_sim=근로소득*(1+1.340+0.019) if prob>=0.9746 & prob<0.9889
	replace 근로_sim=근로소득*(1+1.586+0.019) if prob>=0.9889
	
	replace 사업_sim=(사업소득+재산소득)*(1-0.613+0.080) if prob<0.0204
	replace 사업_sim=(사업소득+재산소득)*(1-0.573+0.080) if prob>=0.0204 & prob<0.0442
	replace 사업_sim=(사업소득+재산소득)*(1-0.528+0.080) if prob>=0.0442 & prob<0.0649
	replace 사업_sim=(사업소득+재산소득)*(1-0.478+0.080) if prob>=0.0649 & prob<0.1043
	replace 사업_sim=(사업소득+재산소득)*(1-0.423+0.080) if prob>=0.1043 & prob<0.1371
	replace 사업_sim=(사업소득+재산소득)*(1-0.362+0.080) if prob>=0.1371 & prob<0.1778
	replace 사업_sim=(사업소득+재산소득)*(1-0.295+0.080) if prob>=0.1778 & prob<0.2242
	replace 사업_sim=(사업소득+재산소득)*(1-0.221+0.080) if prob>=0.2242 & prob<0.2823
	replace 사업_sim=(사업소득+재산소득)*(1-0.139+0.080) if prob>=0.2823 & prob<0.3520
	replace 사업_sim=(사업소득+재산소득)*(1-0.049+0.080) if prob>=0.3520 & prob<0.4292
	replace 사업_sim=(사업소득+재산소득)*(1+0.051+0.080) if prob>=0.4292 & prob<0.5646
	replace 사업_sim=(사업소득+재산소득)*(1+0.162+0.080) if prob>=0.5646 & prob<0.6434
	replace 사업_sim=(사업소득+재산소득)*(1+0.284+0.080) if prob>=0.6434 & prob<0.7179
	replace 사업_sim=(사업소득+재산소득)*(1+0.419+0.080) if prob>=0.7179 & prob<0.7730
	replace 사업_sim=(사업소득+재산소득)*(1+0.568+0.080) if prob>=0.7730 & prob<0.8330
	replace 사업_sim=(사업소득+재산소득)*(1+0.733+0.080) if prob>=0.8330 & prob<0.8777
	replace 사업_sim=(사업소득+재산소득)*(1+0.916+0.080) if prob>=0.8777 & prob<0.9219
	replace 사업_sim=(사업소득+재산소득)*(1+1.117+0.080) if prob>=0.9219 & prob<0.9471
	replace 사업_sim=(사업소득+재산소득)*(1+1.340+0.080) if prob>=0.9471 & prob<0.9730
	replace 사업_sim=(사업소득+재산소득)*(1+1.586+0.080) if prob>=0.9730
	
	replace FM_sim = 근로_sim+사업_sim+이전소득-비소비지출 - DP_sim - BLC_sim
	replace PD_sim = 0 if FM_sim>=0
	replace PD_sim = 0 if (abs(FM_sim)/12)*0.69 <= LA_sim & FM_sim<0
	replace PD_sim = 1- LA_sim/((abs(FM_sim)/12)*0.69) if (abs(FM_sim)/12)*0.69 > LA_sim & FM_sim<0
	replace LGD_ratio_sim = PD_sim*max((부채-CA*0.8)/부채, 0)
	
	summ PD_sim [aweight=가중값]
	replace total_PD = r(mean) in `i'
	summ PD_sim [aweight=부채가중]
	replace total_EAD = r(mean) in `i'
	summ LGD_ratio_sim [aweight=부채가중]
	replace total_LGD = r(mean) in `i'
	summ PD_sim [aweight=가중값] if 은행대출>0
	replace bank_PD = r(mean) in `i'
	summ PD_sim [aweight=은행가중]
	replace bank_EAD = r(mean) in `i'
	summ LGD_ratio_sim [aweight=은행가중]
	replace bank_LGD = r(mean) in `i'
	summ PD_sim [aweight=가중값] if 비은행대출>0
	replace non_PD = r(mean) in `i'
	summ PD_sim [aweight=비은행가중]
	replace non_EAD = r(mean) in `i'
	summ LGD_ratio_sim [aweight=비은행가중]
	replace non_LGD = r(mean) in `i'
}


***** base 시나리오 분석

** base 금리충격 396bp
gen 이자부담증가분_base=최초대출금액*float_exposure*0.0396
gen DP_base=DP+이자부담증가분_base

** base 물가충격 3.7%
gen BLC_base=736.219 if 가구원수==1
replace BLC_base=1237.575 if 가구원수==2
replace BLC_base=1597.379 if 가구원수==3
replace BLC_base=1941.685 if 가구원수==4
replace BLC_base=2227.204 if 가구원수==5
replace BLC_base=2634.395 if 가구원수==6
replace BLC_base=2991.166 if 가구원수==7
replace BLC_base=2991.166+(336.830*(가구원수-7)) if 가구원수>7

** base 주가충격/부동산충격
gen 주식채권펀드_base=주식채권펀드*(1-0.213)*(1+0.121)
gen LA_base=(LA-주식채권펀드)+주식채권펀드_base
gen CA_base=부동산*(1.044)*(1.028)

** base 소득충격 (피보 4.0%)
gen total_PD_base=.
gen total_EAD_base=.
gen total_LGD_base=.
gen bank_PD_base=.
gen bank_EAD_base=.
gen bank_LGD_base=.
gen non_PD_base=.
gen non_EAD_base=.
gen non_LGD_base=.


forvalues i=1(1)1000 {
	set seed `i*7'
	replace prob=runiform()
	
	replace 근로_sim=근로소득*(1-0.713+0.014) if prob<0.0089
	replace 근로_sim=근로소득*(1-0.683+0.014) if prob>=0.0089 & prob<0.0183
	replace 근로_sim=근로소득*(1-0.650+0.014) if prob>=0.0183 & prob<0.0289
	replace 근로_sim=근로소득*(1-0.613+0.014) if prob>=0.0289 & prob<0.0440
	replace 근로_sim=근로소득*(1-0.573+0.014) if prob>=0.0440 & prob<0.0580
	replace 근로_sim=근로소득*(1-0.528+0.014) if prob>=0.0580 & prob<0.0788
	replace 근로_sim=근로소득*(1-0.478+0.014) if prob>=0.0788 & prob<0.0949
	replace 근로_sim=근로소득*(1-0.423+0.014) if prob>=0.0949 & prob<0.1162
	replace 근로_sim=근로소득*(1-0.362+0.014) if prob>=0.1162 & prob<0.1454	
	replace 근로_sim=근로소득*(1-0.295+0.014) if prob>=0.1454 & prob<0.1874	
	replace 근로_sim=근로소득*(1-0.221+0.014) if prob>=0.1874 & prob<0.2247	
	replace 근로_sim=근로소득*(1-0.139+0.014) if prob>=0.2247 & prob<0.2802	
	replace 근로_sim=근로소득*(1-0.049+0.014) if prob>=0.2802 & prob<0.3552	
	replace 근로_sim=근로소득*(1+0.051+0.014) if prob>=0.3552 & prob<0.4706	
	replace 근로_sim=근로소득*(1+0.162+0.014) if prob>=0.4706 & prob<0.6122	
	replace 근로_sim=근로소득*(1+0.284+0.014) if prob>=0.6122 & prob<0.7090	
	replace 근로_sim=근로소득*(1+0.419+0.014) if prob>=0.7090 & prob<0.7767	
	replace 근로_sim=근로소득*(1+0.568+0.014) if prob>=0.7767 & prob<0.8293	
	replace 근로_sim=근로소득*(1+0.733+0.014) if prob>=0.8293 & prob<0.8652	
	replace 근로_sim=근로소득*(1+0.916+0.014) if prob>=0.8652 & prob<0.8997	
	replace 근로_sim=근로소득*(1+1.117+0.014) if prob>=0.8997 & prob<0.9240	
	replace 근로_sim=근로소득*(1+1.340+0.014) if prob>=0.9240 & prob<0.9442	
	replace 근로_sim=근로소득*(1+1.586+0.014) if prob>=0.9442 & prob<0.9614	
	replace 근로_sim=근로소득*(1+1.858+0.014) if prob>=0.9614 & prob<0.9755	
	replace 근로_sim=근로소득*(1+2.158+0.014) if prob>=0.9755 & prob<0.9884	
	replace 근로_sim=근로소득*(1+2.490+0.014) if prob>=0.9884 
	
	replace 사업_sim=사업소득*(1-0.713+0.104) if prob<0.0146
	replace 사업_sim=사업소득*(1-0.683+0.104) if prob>=0.0146 & prob<0.0268
	replace 사업_sim=사업소득*(1-0.650+0.104) if prob>=0.0268 & prob<0.0471
	replace 사업_sim=사업소득*(1-0.613+0.104) if prob>=0.0471 & prob<0.0679
	replace 사업_sim=사업소득*(1-0.573+0.104) if prob>=0.0679 & prob<0.0949
	replace 사업_sim=사업소득*(1-0.528+0.104) if prob>=0.0949 & prob<0.1126
	replace 사업_sim=사업소득*(1-0.478+0.104) if prob>=0.1126 & prob<0.1469
	replace 사업_sim=사업소득*(1-0.423+0.104) if prob>=0.1469 & prob<0.1856
	replace 사업_sim=사업소득*(1-0.362+0.104) if prob>=0.1856 & prob<0.2215	
	replace 사업_sim=사업소득*(1-0.295+0.104) if prob>=0.2215 & prob<0.2581	
	replace 사업_sim=사업소득*(1-0.221+0.104) if prob>=0.2581 & prob<0.3047	
	replace 사업_sim=사업소득*(1-0.139+0.104) if prob>=0.3047 & prob<0.3663	
	replace 사업_sim=사업소득*(1-0.049+0.104) if prob>=0.3663 & prob<0.4305	
	replace 사업_sim=사업소득*(1+0.051+0.104) if prob>=0.4305 & prob<0.5407	
	replace 사업_sim=사업소득*(1+0.162+0.104) if prob>=0.5407 & prob<0.6129	
	replace 사업_sim=사업소득*(1+0.284+0.104) if prob>=0.6129 & prob<0.6779	
	replace 사업_sim=사업소득*(1+0.419+0.104) if prob>=0.6779 & prob<0.7281	
	replace 사업_sim=사업소득*(1+0.568+0.104) if prob>=0.7281 & prob<0.7856	
	replace 사업_sim=사업소득*(1+0.733+0.104) if prob>=0.7856 & prob<0.8245	
	replace 사업_sim=사업소득*(1+0.916+0.104) if prob>=0.8245 & prob<0.8667	
	replace 사업_sim=사업소득*(1+1.117+0.104) if prob>=0.8667 & prob<0.8942	
	replace 사업_sim=사업소득*(1+1.340+0.104) if prob>=0.8942 & prob<0.9179	
	replace 사업_sim=사업소득*(1+1.586+0.104) if prob>=0.9179 & prob<0.9454	
	replace 사업_sim=사업소득*(1+1.858+0.104) if prob>=0.9454 & prob<0.9691	
	replace 사업_sim=사업소득*(1+2.158+0.104) if prob>=0.9691 & prob<0.9847	
	replace 사업_sim=사업소득*(1+2.490+0.104) if prob>=0.9847
	
	replace FM_sim = 근로_sim+사업_sim+이전소득-비소비지출 - DP_base - BLC_base
	replace PD_sim = 0 if FM_sim>=0
	replace PD_sim = 0 if (abs(FM_sim)/12)*0.69 <= LA_base & FM_sim<0
	replace PD_sim = 1- LA_base/((abs(FM_sim)/12)*0.69) if (abs(FM_sim)/12)*0.69 > LA_base & FM_sim<0
	replace LGD_ratio_sim = PD_sim*max((부채-CA_base*0.8)/부채, 0)
	
	summ PD_sim [aweight=가중값]
	replace total_PD_base = r(mean) in `i'
	summ PD_sim [aweight=부채가중]
	replace total_EAD_base = r(mean) in `i'
	summ LGD_ratio_sim [aweight=부채가중]
	replace total_LGD_base = r(mean) in `i'
	summ PD_sim [aweight=가중값] if 은행대출>0
	replace bank_PD_base = r(mean) in `i'
	summ PD_sim [aweight=은행가중]
	replace bank_EAD_base = r(mean) in `i'
	summ LGD_ratio_sim [aweight=은행가중]
	replace bank_LGD_base = r(mean) in `i'
	summ PD_sim [aweight=가중값] if 비은행대출>0
	replace non_PD_base = r(mean) in `i'
	summ PD_sim [aweight=비은행가중]
	replace non_EAD_base = r(mean) in `i'
	summ LGD_ratio_sim [aweight=비은행가중]
	replace non_LGD_base = r(mean) in `i'
}




***** weak 시나리오 분석

** weak 금리충격 321bp
gen 이자부담증가분_weak=최초대출금액*float_exposure*0.0321
gen DP_weak=DP+이자부담증가분_weak

** weak 물가충격 2.7%
gen BLC_weak=731.862 if 가구원수==1
replace BLC_weak=1231.025 if 가구원수==2
replace BLC_weak=1588.733 if 가구원수==3
replace BLC_weak=1933.973 if 가구원수==4
replace BLC_weak=2228.478 if 가구원수==5
replace BLC_weak=2623.188 if 가구원수==6
replace BLC_weak=2971.075 if 가구원수==7
replace BLC_weak=2971.075+(336.547*(가구원수-7)) if 가구원수>7

** weak 주가충격/부동산충격
gen 주식채권펀드_weak=주식채권펀드*(1-0.213)*(1+0.382)
gen LA_weak=(LA-주식채권펀드)+주식채권펀드_weak
gen CA_weak=부동산*(1.044)*(1.061)


** weak 소득충격 
gen total_PD_weak=.
gen total_EAD_weak=.
gen total_LGD_weak=.
gen bank_PD_weak=.
gen bank_EAD_weak=.
gen bank_LGD_weak=.
gen non_PD_weak=.
gen non_EAD_weak=.
gen non_LGD_weak=.


forvalues i=1(1)1000 {
	set seed `i*7'
	replace prob=runiform()
	
	replace 근로_sim=근로소득*(1-0.713+0.036) if prob<0.0089
	replace 근로_sim=근로소득*(1-0.683+0.036) if prob>=0.0089 & prob<0.0183
	replace 근로_sim=근로소득*(1-0.650+0.036) if prob>=0.0183 & prob<0.0289
	replace 근로_sim=근로소득*(1-0.613+0.036) if prob>=0.0289 & prob<0.0440
	replace 근로_sim=근로소득*(1-0.573+0.036) if prob>=0.0440 & prob<0.0580
	replace 근로_sim=근로소득*(1-0.528+0.036) if prob>=0.0580 & prob<0.0788
	replace 근로_sim=근로소득*(1-0.478+0.036) if prob>=0.0788 & prob<0.0949
	replace 근로_sim=근로소득*(1-0.423+0.036) if prob>=0.0949 & prob<0.1162
	replace 근로_sim=근로소득*(1-0.362+0.036) if prob>=0.1162 & prob<0.1454	
	replace 근로_sim=근로소득*(1-0.295+0.036) if prob>=0.1454 & prob<0.1874	
	replace 근로_sim=근로소득*(1-0.221+0.036) if prob>=0.1874 & prob<0.2247	
	replace 근로_sim=근로소득*(1-0.139+0.036) if prob>=0.2247 & prob<0.2802	
	replace 근로_sim=근로소득*(1-0.049+0.036) if prob>=0.2802 & prob<0.3552	
	replace 근로_sim=근로소득*(1+0.051+0.036) if prob>=0.3552 & prob<0.4706	
	replace 근로_sim=근로소득*(1+0.162+0.036) if prob>=0.4706 & prob<0.6122	
	replace 근로_sim=근로소득*(1+0.284+0.036) if prob>=0.6122 & prob<0.7090	
	replace 근로_sim=근로소득*(1+0.419+0.036) if prob>=0.7090 & prob<0.7767	
	replace 근로_sim=근로소득*(1+0.568+0.036) if prob>=0.7767 & prob<0.8293	
	replace 근로_sim=근로소득*(1+0.733+0.036) if prob>=0.8293 & prob<0.8652	
	replace 근로_sim=근로소득*(1+0.916+0.036) if prob>=0.8652 & prob<0.8997	
	replace 근로_sim=근로소득*(1+1.117+0.036) if prob>=0.8997 & prob<0.9240	
	replace 근로_sim=근로소득*(1+1.340+0.036) if prob>=0.9240 & prob<0.9442	
	replace 근로_sim=근로소득*(1+1.586+0.036) if prob>=0.9442 & prob<0.9614	
	replace 근로_sim=근로소득*(1+1.858+0.036) if prob>=0.9614 & prob<0.9755	
	replace 근로_sim=근로소득*(1+2.158+0.036) if prob>=0.9755 & prob<0.9884	
	replace 근로_sim=근로소득*(1+2.490+0.036) if prob>=0.9884 
	
	replace 사업_sim=사업소득*(1-0.713+0.126) if prob<0.0146
	replace 사업_sim=사업소득*(1-0.683+0.126) if prob>=0.0146 & prob<0.0268
	replace 사업_sim=사업소득*(1-0.650+0.126) if prob>=0.0268 & prob<0.0471
	replace 사업_sim=사업소득*(1-0.613+0.126) if prob>=0.0471 & prob<0.0679
	replace 사업_sim=사업소득*(1-0.573+0.126) if prob>=0.0679 & prob<0.0949
	replace 사업_sim=사업소득*(1-0.528+0.126) if prob>=0.0949 & prob<0.1126
	replace 사업_sim=사업소득*(1-0.478+0.126) if prob>=0.1126 & prob<0.1469
	replace 사업_sim=사업소득*(1-0.423+0.126) if prob>=0.1469 & prob<0.1856
	replace 사업_sim=사업소득*(1-0.362+0.126) if prob>=0.1856 & prob<0.2215	
	replace 사업_sim=사업소득*(1-0.295+0.126) if prob>=0.2215 & prob<0.2581	
	replace 사업_sim=사업소득*(1-0.221+0.126) if prob>=0.2581 & prob<0.3047	
	replace 사업_sim=사업소득*(1-0.139+0.126) if prob>=0.3047 & prob<0.3663	
	replace 사업_sim=사업소득*(1-0.049+0.126) if prob>=0.3663 & prob<0.4305	
	replace 사업_sim=사업소득*(1+0.051+0.126) if prob>=0.4305 & prob<0.5407	
	replace 사업_sim=사업소득*(1+0.162+0.126) if prob>=0.5407 & prob<0.6129	
	replace 사업_sim=사업소득*(1+0.284+0.126) if prob>=0.6129 & prob<0.6779	
	replace 사업_sim=사업소득*(1+0.419+0.126) if prob>=0.6779 & prob<0.7281	
	replace 사업_sim=사업소득*(1+0.568+0.126) if prob>=0.7281 & prob<0.7856	
	replace 사업_sim=사업소득*(1+0.733+0.126) if prob>=0.7856 & prob<0.8245	
	replace 사업_sim=사업소득*(1+0.916+0.126) if prob>=0.8245 & prob<0.8667	
	replace 사업_sim=사업소득*(1+1.117+0.126) if prob>=0.8667 & prob<0.8942	
	replace 사업_sim=사업소득*(1+1.340+0.126) if prob>=0.8942 & prob<0.9179	
	replace 사업_sim=사업소득*(1+1.586+0.126) if prob>=0.9179 & prob<0.9454	
	replace 사업_sim=사업소득*(1+1.858+0.126) if prob>=0.9454 & prob<0.9691	
	replace 사업_sim=사업소득*(1+2.158+0.126) if prob>=0.9691 & prob<0.9847	
	replace 사업_sim=사업소득*(1+2.490+0.126) if prob>=0.9847
	
	replace FM_sim = 근로_sim+사업_sim+이전소득-비소비지출 - DP_weak - BLC_weak
	replace PD_sim = 0 if FM_sim>=0
	replace PD_sim = 0 if (abs(FM_sim)/12)*0.69 <= LA_weak & FM_sim<0
	replace PD_sim = 1- LA_weak/((abs(FM_sim)/12)*0.69) if (abs(FM_sim)/12)*0.69 > LA_weak & FM_sim<0
	replace LGD_ratio_sim = PD_sim*max((부채-CA_weak*0.8)/부채, 0)
	
	summ PD_sim [aweight=가중값]
	replace total_PD_weak = r(mean) in `i'
	summ PD_sim [aweight=부채가중]
	replace total_EAD_weak = r(mean) in `i'
	summ LGD_ratio_sim [aweight=부채가중]
	replace total_LGD_weak = r(mean) in `i'
	summ PD_sim [aweight=가중값] if 은행대출>0
	replace bank_PD_weak = r(mean) in `i'
	summ PD_sim [aweight=은행가중]
	replace bank_EAD_weak = r(mean) in `i'
	summ LGD_ratio_sim [aweight=은행가중]
	replace bank_LGD_weak = r(mean) in `i'
	summ PD_sim [aweight=가중값] if 비은행대출>0
	replace non_PD_weak = r(mean) in `i'
	summ PD_sim [aweight=비은행가중]
	replace non_EAD_weak = r(mean) in `i'
	summ LGD_ratio_sim [aweight=비은행가중]
	replace non_LGD_weak = r(mean) in `i'
}








***** extreme 시나리오 분석

** extreme 금리충격 454bp
gen 이자부담증가분_extreme=최초대출금액*float_exposure*0.0454
gen DP_extreme=DP+이자부담증가분_extreme

** extreme 물가충격 6%
gen BLC_extreme=740.577 if 가구원수==1
replace BLC_extreme=1244.125 if 가구원수==2
replace BLC_extreme=1606.025 if 가구원수==3
replace BLC_extreme=1949.398 if 가구원수==4
replace BLC_extreme=2225.929 if 가구원수==5
replace BLC_extreme=2645.601 if 가구원수==6
replace BLC_extreme=3011.258 if 가구원수==7
replace BLC_extreme=3011.258+(337.112*(가구원수-7)) if 가구원수>7

** extreme 주가충격/부동산충격
gen 주식채권펀드_extreme=주식채권펀드*(1-0.213)*(1-0.141)
gen LA_extreme=(LA-주식채권펀드)+주식채권펀드_extreme
gen CA_extreme=부동산*(1.044)*(1-0.005)


** extreme 소득충격 
gen total_PD_extreme=.
gen total_EAD_extreme=.
gen total_LGD_extreme=.
gen bank_PD_extreme=.
gen bank_EAD_extreme=.
gen bank_LGD_extreme=.
gen non_PD_extreme=.
gen non_EAD_extreme=.
gen non_LGD_extreme=.


forvalues i=1(1)1000 {
	set seed `i*7'
	replace prob=runiform()
	
	replace 근로_sim=근로소득*(1-0.713-0.008) if prob<0.0089
	replace 근로_sim=근로소득*(1-0.683-0.008) if prob>=0.0089 & prob<0.0183
	replace 근로_sim=근로소득*(1-0.650-0.008) if prob>=0.0183 & prob<0.0289
	replace 근로_sim=근로소득*(1-0.613-0.008) if prob>=0.0289 & prob<0.0440
	replace 근로_sim=근로소득*(1-0.573-0.008) if prob>=0.0440 & prob<0.0580
	replace 근로_sim=근로소득*(1-0.528-0.008) if prob>=0.0580 & prob<0.0788
	replace 근로_sim=근로소득*(1-0.478-0.008) if prob>=0.0788 & prob<0.0949
	replace 근로_sim=근로소득*(1-0.423-0.008) if prob>=0.0949 & prob<0.1162
	replace 근로_sim=근로소득*(1-0.362-0.008) if prob>=0.1162 & prob<0.1454	
	replace 근로_sim=근로소득*(1-0.295-0.008) if prob>=0.1454 & prob<0.1874	
	replace 근로_sim=근로소득*(1-0.221-0.008) if prob>=0.1874 & prob<0.2247	
	replace 근로_sim=근로소득*(1-0.139-0.008) if prob>=0.2247 & prob<0.2802	
	replace 근로_sim=근로소득*(1-0.049-0.008) if prob>=0.2802 & prob<0.3552	
	replace 근로_sim=근로소득*(1+0.051-0.008) if prob>=0.3552 & prob<0.4706	
	replace 근로_sim=근로소득*(1+0.162-0.008) if prob>=0.4706 & prob<0.6122	
	replace 근로_sim=근로소득*(1+0.284-0.008) if prob>=0.6122 & prob<0.7090	
	replace 근로_sim=근로소득*(1+0.419-0.008) if prob>=0.7090 & prob<0.7767	
	replace 근로_sim=근로소득*(1+0.568-0.008) if prob>=0.7767 & prob<0.8293	
	replace 근로_sim=근로소득*(1+0.733-0.008) if prob>=0.8293 & prob<0.8652	
	replace 근로_sim=근로소득*(1+0.916-0.008) if prob>=0.8652 & prob<0.8997	
	replace 근로_sim=근로소득*(1+1.117-0.008) if prob>=0.8997 & prob<0.9240	
	replace 근로_sim=근로소득*(1+1.340-0.008) if prob>=0.9240 & prob<0.9442	
	replace 근로_sim=근로소득*(1+1.586-0.008) if prob>=0.9442 & prob<0.9614	
	replace 근로_sim=근로소득*(1+1.858-0.008) if prob>=0.9614 & prob<0.9755	
	replace 근로_sim=근로소득*(1+2.158-0.008) if prob>=0.9755 & prob<0.9884	
	replace 근로_sim=근로소득*(1+2.490-0.008) if prob>=0.9884 
	
	replace 사업_sim=사업소득*(1-0.713+0.081) if prob<0.0146
	replace 사업_sim=사업소득*(1-0.683+0.081) if prob>=0.0146 & prob<0.0268
	replace 사업_sim=사업소득*(1-0.650+0.081) if prob>=0.0268 & prob<0.0471
	replace 사업_sim=사업소득*(1-0.613+0.081) if prob>=0.0471 & prob<0.0679
	replace 사업_sim=사업소득*(1-0.573+0.081) if prob>=0.0679 & prob<0.0949
	replace 사업_sim=사업소득*(1-0.528+0.081) if prob>=0.0949 & prob<0.1126
	replace 사업_sim=사업소득*(1-0.478+0.081) if prob>=0.1126 & prob<0.1469
	replace 사업_sim=사업소득*(1-0.423+0.081) if prob>=0.1469 & prob<0.1856
	replace 사업_sim=사업소득*(1-0.362+0.081) if prob>=0.1856 & prob<0.2215	
	replace 사업_sim=사업소득*(1-0.295+0.081) if prob>=0.2215 & prob<0.2581	
	replace 사업_sim=사업소득*(1-0.221+0.081) if prob>=0.2581 & prob<0.3047	
	replace 사업_sim=사업소득*(1-0.139+0.081) if prob>=0.3047 & prob<0.3663	
	replace 사업_sim=사업소득*(1-0.049+0.081) if prob>=0.3663 & prob<0.4305	
	replace 사업_sim=사업소득*(1+0.051+0.081) if prob>=0.4305 & prob<0.5407	
	replace 사업_sim=사업소득*(1+0.162+0.081) if prob>=0.5407 & prob<0.6129	
	replace 사업_sim=사업소득*(1+0.284+0.081) if prob>=0.6129 & prob<0.6779	
	replace 사업_sim=사업소득*(1+0.419+0.081) if prob>=0.6779 & prob<0.7281	
	replace 사업_sim=사업소득*(1+0.568+0.081) if prob>=0.7281 & prob<0.7856	
	replace 사업_sim=사업소득*(1+0.733+0.081) if prob>=0.7856 & prob<0.8245	
	replace 사업_sim=사업소득*(1+0.916+0.081) if prob>=0.8245 & prob<0.8667	
	replace 사업_sim=사업소득*(1+1.117+0.081) if prob>=0.8667 & prob<0.8942	
	replace 사업_sim=사업소득*(1+1.340+0.081) if prob>=0.8942 & prob<0.9179	
	replace 사업_sim=사업소득*(1+1.586+0.081) if prob>=0.9179 & prob<0.9454	
	replace 사업_sim=사업소득*(1+1.858+0.081) if prob>=0.9454 & prob<0.9691	
	replace 사업_sim=사업소득*(1+2.158+0.081) if prob>=0.9691 & prob<0.9847	
	replace 사업_sim=사업소득*(1+2.490+0.081) if prob>=0.9847
	
	replace FM_sim = 근로_sim+사업_sim+이전소득-비소비지출 - DP_extreme - BLC_extreme
	replace PD_sim = 0 if FM_sim>=0
	replace PD_sim = 0 if (abs(FM_sim)/12)*0.69 <= LA_extreme & FM_sim<0
	replace PD_sim = 1- LA_extreme/((abs(FM_sim)/12)*0.69) if (abs(FM_sim)/12)*0.69 > LA_extreme & FM_sim<0
	replace LGD_ratio_sim = PD_sim*max((부채-CA_extreme*0.8)/부채, 0)
	
	summ PD_sim [aweight=가중값]
	replace total_PD_extreme = r(mean) in `i'
	summ PD_sim [aweight=부채가중]
	replace total_EAD_extreme = r(mean) in `i'
	summ LGD_ratio_sim [aweight=부채가중]
	replace total_LGD_extreme = r(mean) in `i'
	summ PD_sim [aweight=가중값] if 은행대출>0
	replace bank_PD_extreme = r(mean) in `i'
	summ PD_sim [aweight=은행가중]
	replace bank_EAD_extreme = r(mean) in `i'
	summ LGD_ratio_sim [aweight=은행가중]
	replace bank_LGD_extreme = r(mean) in `i'
	summ PD_sim [aweight=가중값] if 비은행대출>0
	replace non_PD_extreme = r(mean) in `i'
	summ PD_sim [aweight=비은행가중]
	replace non_EAD_extreme = r(mean) in `i'
	summ LGD_ratio_sim [aweight=비은행가중]
	replace non_LGD_extreme = r(mean) in `i'
}




***** extension

** base 금리충격 0부터 1000bp
replace 이자부담증가분_base=최초대출금액*float_exposure*0.10
replace DP_base=DP+이자부담증가분_base


forvalues i=1(1)1000 {
	set seed `i*7'
	replace prob=runiform()
	
	replace 근로_sim=근로소득*(1-0.713+0.014) if prob<0.0089
	replace 근로_sim=근로소득*(1-0.683+0.014) if prob>=0.0089 & prob<0.0183
	replace 근로_sim=근로소득*(1-0.650+0.014) if prob>=0.0183 & prob<0.0289
	replace 근로_sim=근로소득*(1-0.613+0.014) if prob>=0.0289 & prob<0.0440
	replace 근로_sim=근로소득*(1-0.573+0.014) if prob>=0.0440 & prob<0.0580
	replace 근로_sim=근로소득*(1-0.528+0.014) if prob>=0.0580 & prob<0.0788
	replace 근로_sim=근로소득*(1-0.478+0.014) if prob>=0.0788 & prob<0.0949
	replace 근로_sim=근로소득*(1-0.423+0.014) if prob>=0.0949 & prob<0.1162
	replace 근로_sim=근로소득*(1-0.362+0.014) if prob>=0.1162 & prob<0.1454	
	replace 근로_sim=근로소득*(1-0.295+0.014) if prob>=0.1454 & prob<0.1874	
	replace 근로_sim=근로소득*(1-0.221+0.014) if prob>=0.1874 & prob<0.2247	
	replace 근로_sim=근로소득*(1-0.139+0.014) if prob>=0.2247 & prob<0.2802	
	replace 근로_sim=근로소득*(1-0.049+0.014) if prob>=0.2802 & prob<0.3552	
	replace 근로_sim=근로소득*(1+0.051+0.014) if prob>=0.3552 & prob<0.4706	
	replace 근로_sim=근로소득*(1+0.162+0.014) if prob>=0.4706 & prob<0.6122	
	replace 근로_sim=근로소득*(1+0.284+0.014) if prob>=0.6122 & prob<0.7090	
	replace 근로_sim=근로소득*(1+0.419+0.014) if prob>=0.7090 & prob<0.7767	
	replace 근로_sim=근로소득*(1+0.568+0.014) if prob>=0.7767 & prob<0.8293	
	replace 근로_sim=근로소득*(1+0.733+0.014) if prob>=0.8293 & prob<0.8652	
	replace 근로_sim=근로소득*(1+0.916+0.014) if prob>=0.8652 & prob<0.8997	
	replace 근로_sim=근로소득*(1+1.117+0.014) if prob>=0.8997 & prob<0.9240	
	replace 근로_sim=근로소득*(1+1.340+0.014) if prob>=0.9240 & prob<0.9442	
	replace 근로_sim=근로소득*(1+1.586+0.014) if prob>=0.9442 & prob<0.9614	
	replace 근로_sim=근로소득*(1+1.858+0.014) if prob>=0.9614 & prob<0.9755	
	replace 근로_sim=근로소득*(1+2.158+0.014) if prob>=0.9755 & prob<0.9884	
	replace 근로_sim=근로소득*(1+2.490+0.014) if prob>=0.9884 
	
	replace 사업_sim=사업소득*(1-0.713+0.104) if prob<0.0146
	replace 사업_sim=사업소득*(1-0.683+0.104) if prob>=0.0146 & prob<0.0268
	replace 사업_sim=사업소득*(1-0.650+0.104) if prob>=0.0268 & prob<0.0471
	replace 사업_sim=사업소득*(1-0.613+0.104) if prob>=0.0471 & prob<0.0679
	replace 사업_sim=사업소득*(1-0.573+0.104) if prob>=0.0679 & prob<0.0949
	replace 사업_sim=사업소득*(1-0.528+0.104) if prob>=0.0949 & prob<0.1126
	replace 사업_sim=사업소득*(1-0.478+0.104) if prob>=0.1126 & prob<0.1469
	replace 사업_sim=사업소득*(1-0.423+0.104) if prob>=0.1469 & prob<0.1856
	replace 사업_sim=사업소득*(1-0.362+0.104) if prob>=0.1856 & prob<0.2215	
	replace 사업_sim=사업소득*(1-0.295+0.104) if prob>=0.2215 & prob<0.2581	
	replace 사업_sim=사업소득*(1-0.221+0.104) if prob>=0.2581 & prob<0.3047	
	replace 사업_sim=사업소득*(1-0.139+0.104) if prob>=0.3047 & prob<0.3663	
	replace 사업_sim=사업소득*(1-0.049+0.104) if prob>=0.3663 & prob<0.4305	
	replace 사업_sim=사업소득*(1+0.051+0.104) if prob>=0.4305 & prob<0.5407	
	replace 사업_sim=사업소득*(1+0.162+0.104) if prob>=0.5407 & prob<0.6129	
	replace 사업_sim=사업소득*(1+0.284+0.104) if prob>=0.6129 & prob<0.6779	
	replace 사업_sim=사업소득*(1+0.419+0.104) if prob>=0.6779 & prob<0.7281	
	replace 사업_sim=사업소득*(1+0.568+0.104) if prob>=0.7281 & prob<0.7856	
	replace 사업_sim=사업소득*(1+0.733+0.104) if prob>=0.7856 & prob<0.8245	
	replace 사업_sim=사업소득*(1+0.916+0.104) if prob>=0.8245 & prob<0.8667	
	replace 사업_sim=사업소득*(1+1.117+0.104) if prob>=0.8667 & prob<0.8942	
	replace 사업_sim=사업소득*(1+1.340+0.104) if prob>=0.8942 & prob<0.9179	
	replace 사업_sim=사업소득*(1+1.586+0.104) if prob>=0.9179 & prob<0.9454	
	replace 사업_sim=사업소득*(1+1.858+0.104) if prob>=0.9454 & prob<0.9691	
	replace 사업_sim=사업소득*(1+2.158+0.104) if prob>=0.9691 & prob<0.9847	
	replace 사업_sim=사업소득*(1+2.490+0.104) if prob>=0.9847
	
	replace FM_sim = 근로_sim+사업_sim+이전소득-비소비지출 - DP_base - BLC_base
	replace PD_sim = 0 if FM_sim>=0
	replace PD_sim = 0 if (abs(FM_sim)/12)*0.69 <= LA_base & FM_sim<0
	replace PD_sim = 1- LA_base/((abs(FM_sim)/12)*0.69) if (abs(FM_sim)/12)*0.69 > LA_base & FM_sim<0
	replace LGD_ratio_sim = PD_sim*max((부채-CA_base*0.8)/부채, 0)
	
	summ PD_sim [aweight=가중값]
	replace total_PD_base = r(mean) in `i'
	summ PD_sim [aweight=부채가중]
	replace total_EAD_base = r(mean) in `i'
	summ LGD_ratio_sim [aweight=부채가중]
	replace total_LGD_base = r(mean) in `i'
	summ PD_sim [aweight=가중값] if 은행대출>0
	replace bank_PD_base = r(mean) in `i'
	summ PD_sim [aweight=은행가중]
	replace bank_EAD_base = r(mean) in `i'
	summ LGD_ratio_sim [aweight=은행가중]
	replace bank_LGD_base = r(mean) in `i'
	summ PD_sim [aweight=가중값] if 비은행대출>0
	replace non_PD_base = r(mean) in `i'
	summ PD_sim [aweight=비은행가중]
	replace non_EAD_base = r(mean) in `i'
	summ LGD_ratio_sim [aweight=비은행가중]
	replace non_LGD_base = r(mean) in `i'
}


****** 소득충격 변화

** base 금리충격 396bp
replace 이자부담증가분_base=최초대출금액*float_exposure*0.0396
replace DP_base=DP+이자부담증가분_base

** base 물가충격 변화
replace BLC_base=736.219 if 가구원수==1
replace BLC_base=1237.575 if 가구원수==2
replace BLC_base=1597.379 if 가구원수==3
replace BLC_base=1941.685 if 가구원수==4
replace BLC_base=2227.204 if 가구원수==5
replace BLC_base=2634.395 if 가구원수==6
replace BLC_base=2991.166 if 가구원수==7
replace BLC_base=2991.166+(336.830*(가구원수-7)) if 가구원수>7

** base 주가충격/부동산충격
replace 주식채권펀드_base=주식채권펀드*(1-0.213)*(1+0.121)
replace LA_base=(LA-주식채권펀드)+주식채권펀드_base
replace CA_base=부동산*(1.044)*(1+0.127)

** base 소득충격 * 변화

forvalues i=1(1)1000 {
	set seed `i*7'
	replace prob=runiform()
	
	replace 근로_sim=근로소득*(1-0.713+0.014) if prob<0.0089
	replace 근로_sim=근로소득*(1-0.683+0.014) if prob>=0.0089 & prob<0.0183
	replace 근로_sim=근로소득*(1-0.650+0.014) if prob>=0.0183 & prob<0.0289
	replace 근로_sim=근로소득*(1-0.613+0.014) if prob>=0.0289 & prob<0.0440
	replace 근로_sim=근로소득*(1-0.573+0.014) if prob>=0.0440 & prob<0.0580
	replace 근로_sim=근로소득*(1-0.528+0.014) if prob>=0.0580 & prob<0.0788
	replace 근로_sim=근로소득*(1-0.478+0.014) if prob>=0.0788 & prob<0.0949
	replace 근로_sim=근로소득*(1-0.423+0.014) if prob>=0.0949 & prob<0.1162
	replace 근로_sim=근로소득*(1-0.362+0.014) if prob>=0.1162 & prob<0.1454	
	replace 근로_sim=근로소득*(1-0.295+0.014) if prob>=0.1454 & prob<0.1874	
	replace 근로_sim=근로소득*(1-0.221+0.014) if prob>=0.1874 & prob<0.2247	
	replace 근로_sim=근로소득*(1-0.139+0.014) if prob>=0.2247 & prob<0.2802	
	replace 근로_sim=근로소득*(1-0.049+0.014) if prob>=0.2802 & prob<0.3552	
	replace 근로_sim=근로소득*(1+0.051+0.014) if prob>=0.3552 & prob<0.4706	
	replace 근로_sim=근로소득*(1+0.162+0.014) if prob>=0.4706 & prob<0.6122	
	replace 근로_sim=근로소득*(1+0.284+0.014) if prob>=0.6122 & prob<0.7090	
	replace 근로_sim=근로소득*(1+0.419+0.014) if prob>=0.7090 & prob<0.7767	
	replace 근로_sim=근로소득*(1+0.568+0.014) if prob>=0.7767 & prob<0.8293	
	replace 근로_sim=근로소득*(1+0.733+0.014) if prob>=0.8293 & prob<0.8652	
	replace 근로_sim=근로소득*(1+0.916+0.014) if prob>=0.8652 & prob<0.8997	
	replace 근로_sim=근로소득*(1+1.117+0.014) if prob>=0.8997 & prob<0.9240	
	replace 근로_sim=근로소득*(1+1.340+0.014) if prob>=0.9240 & prob<0.9442	
	replace 근로_sim=근로소득*(1+1.586+0.014) if prob>=0.9442 & prob<0.9614	
	replace 근로_sim=근로소득*(1+1.858+0.014) if prob>=0.9614 & prob<0.9755	
	replace 근로_sim=근로소득*(1+2.158+0.014) if prob>=0.9755 & prob<0.9884	
	replace 근로_sim=근로소득*(1+2.490+0.014) if prob>=0.9884 
	
	replace 사업_sim=사업소득*(1-0.713+0.104) if prob<0.0146
	replace 사업_sim=사업소득*(1-0.683+0.104) if prob>=0.0146 & prob<0.0268
	replace 사업_sim=사업소득*(1-0.650+0.104) if prob>=0.0268 & prob<0.0471
	replace 사업_sim=사업소득*(1-0.613+0.104) if prob>=0.0471 & prob<0.0679
	replace 사업_sim=사업소득*(1-0.573+0.104) if prob>=0.0679 & prob<0.0949
	replace 사업_sim=사업소득*(1-0.528+0.104) if prob>=0.0949 & prob<0.1126
	replace 사업_sim=사업소득*(1-0.478+0.104) if prob>=0.1126 & prob<0.1469
	replace 사업_sim=사업소득*(1-0.423+0.104) if prob>=0.1469 & prob<0.1856
	replace 사업_sim=사업소득*(1-0.362+0.104) if prob>=0.1856 & prob<0.2215	
	replace 사업_sim=사업소득*(1-0.295+0.104) if prob>=0.2215 & prob<0.2581	
	replace 사업_sim=사업소득*(1-0.221+0.104) if prob>=0.2581 & prob<0.3047	
	replace 사업_sim=사업소득*(1-0.139+0.104) if prob>=0.3047 & prob<0.3663	
	replace 사업_sim=사업소득*(1-0.049+0.104) if prob>=0.3663 & prob<0.4305	
	replace 사업_sim=사업소득*(1+0.051+0.104) if prob>=0.4305 & prob<0.5407	
	replace 사업_sim=사업소득*(1+0.162+0.104) if prob>=0.5407 & prob<0.6129	
	replace 사업_sim=사업소득*(1+0.284+0.104) if prob>=0.6129 & prob<0.6779	
	replace 사업_sim=사업소득*(1+0.419+0.104) if prob>=0.6779 & prob<0.7281	
	replace 사업_sim=사업소득*(1+0.568+0.104) if prob>=0.7281 & prob<0.7856	
	replace 사업_sim=사업소득*(1+0.733+0.104) if prob>=0.7856 & prob<0.8245	
	replace 사업_sim=사업소득*(1+0.916+0.104) if prob>=0.8245 & prob<0.8667	
	replace 사업_sim=사업소득*(1+1.117+0.104) if prob>=0.8667 & prob<0.8942	
	replace 사업_sim=사업소득*(1+1.340+0.104) if prob>=0.8942 & prob<0.9179	
	replace 사업_sim=사업소득*(1+1.586+0.104) if prob>=0.9179 & prob<0.9454	
	replace 사업_sim=사업소득*(1+1.858+0.104) if prob>=0.9454 & prob<0.9691	
	replace 사업_sim=사업소득*(1+2.158+0.104) if prob>=0.9691 & prob<0.9847	
	replace 사업_sim=사업소득*(1+2.490+0.104) if prob>=0.9847
	
	replace FM_sim = 근로_sim+사업_sim+이전소득-비소비지출 - DP_base - BLC_base
	replace PD_sim = 0 if FM_sim>=0
	replace PD_sim = 0 if (abs(FM_sim)/12)*0.69 <= LA_base & FM_sim<0
	replace PD_sim = 1- LA_base/((abs(FM_sim)/12)*0.69) if (abs(FM_sim)/12)*0.69 > LA_base & FM_sim<0
	replace LGD_ratio_sim = PD_sim*max((부채-CA_base*0.8)/부채, 0)
	
	summ PD_sim [aweight=가중값]
	replace total_PD_base = r(mean) in `i'
	summ PD_sim [aweight=부채가중]
	replace total_EAD_base = r(mean) in `i'
	summ LGD_ratio_sim [aweight=부채가중]
	replace total_LGD_base = r(mean) in `i'
	summ PD_sim [aweight=가중값] if 은행대출>0
	replace bank_PD_base = r(mean) in `i'
	summ PD_sim [aweight=은행가중]
	replace bank_EAD_base = r(mean) in `i'
	summ LGD_ratio_sim [aweight=은행가중]
	replace bank_LGD_base = r(mean) in `i'
	summ PD_sim [aweight=가중값] if 비은행대출>0
	replace non_PD_base = r(mean) in `i'
	summ PD_sim [aweight=비은행가중]
	replace non_EAD_base = r(mean) in `i'
	summ LGD_ratio_sim [aweight=비은행가중]
	replace non_LGD_base = r(mean) in `i'
}






