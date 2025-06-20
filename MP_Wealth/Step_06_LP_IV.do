* Step_06_LP_IV
* Yeongwoong Do (2025)

* This do-file must be run after "Step_01".

clear
cls

cd "C:\Users\82106\Documents\5. 2025\MP_Wealth_inequality\replication_package"

**************<Generate 02_LP_data.dta>*****************************************

* (1) MPS (monthly->quarterly)
import excel "raw_data\monetary-policy-surprises-data.xlsx", sheet("Monthly (update 2023)") firstrow

gen quarter = .
replace quarter=1 if Month<4
replace quarter=2 if Month>=4 & Month<7
replace quarter=3 if Month>=7 & Month<10
replace quarter=4 if Month>=10

collapse (sum) MPS_ORTH, by(Year quarter)
save "STATA_dta\02_LP_data.dta", replace

* (2) Macrovariables

import excel "raw_data\FRED_data.xlsx", sheet("quarter") firstrow clear

merge m:1 Year quarter using STATA_dta\02_LP_data.dta
drop _merge

gen time = tq(1988,1)-1+_n
format time %tq
drop Year quarter
order time

* (3) DFA Gini
merge m:1 time using STATA_dta\01_DFA_gini.dta
drop _merge

rename (DGS1 CPIAUCSL GDPC1 ebp MPS_ORTH ps_gini_DA) (R Pi Y EBP MPS Gini)
save "STATA_dta\02_LP_data.dta", replace

cap mkdir Table
export excel using "Table\T5_dataset.xlsx", firstrow(variables) replace

**************<Run LP-IV>*******************************************************

* (1) Only R and Gini
use STATA_dta\02_LP_data.dta, clear
* cut the sample
drop if time>tq(2023,4)

tsset time
gen dGini = Gini - L.Gini
gen dlnY  = log(Y) - log(L.Y)
gen dlnPi = log(Pi) - log(L.Pi)

local p = 4 

* First stage
reg R MPS L(1/`p').R L(1/`p').dGini
test MPS

* LP-IV
gen beta = .
gen se_w = .
gen se_nw = .
local row = 1
local H   = 12

forvalues h = 0/`H'{
	* dependent variable
	gen Y_dif_`h' = F`h'.Gini - L.Gini
	
	* LP-IV
	qui ivreg2 Y_dif_`h' L(1/`p').R L(1/`p').dGini (R = MPS), r
	* qui ivregress 2sls Y_dif_`h' L(1/`p').R L(1/`p').dGini (R = MPS L(1/`p').R L(1/`p').dGini), r
	
	* results
	replace beta = _b[R] in `row'
	replace se_w   = _se[R] in `row'
	
	* NW
	local NW_h = `h'+1
	qui ivreg2 Y_dif_`h' L(1/`p').R L(1/`p').dGini (R = MPS), bw(`NW_h') r
	replace se_nw   = _se[R] in `row'
		
	local ++row
}

keep beta se_w se_nw

gen horizon=0-1+_n
drop if horizon>`H'
gen LB = beta-se_w
gen UB = beta+se_w
gen LB2 = beta-1.645*se_w
gen UB2 = beta+1.645*se_w

* main
tsset horizon

twoway (rarea UB LB horizon, fcolor(red%30) lcolor(%0)) (rarea UB2 LB2 horizon, fcolor(red%10) lcolor(%0)) (line beta horizon, lcolor(red) lwidth(thick)), yline(0, lcolor(black) lpattern(dash)) xscale(range(0 `H')) xlabel(0(2)`H') graphregion(color(white))  ytitle("Difference of Gini") xtitle("Quarter") legend(off) name(g1, replace)

cap mkdir Figure	
graph export "Figure/F3.png", replace

* NW
gen LB_nw = beta-se_nw
gen UB_nw = beta+se_nw
gen LB2_nw = beta-1.645*se_nw
gen UB2_nw = beta+1.645*se_nw

twoway (rarea UB_nw LB_nw horizon, fcolor(red%30) lcolor(%0)) (rarea UB2_nw LB2_nw horizon, fcolor(red%10) lcolor(%0)) (line beta horizon, lcolor(red) lwidth(thick)), yline(0, lcolor(black) lpattern(dash)) xscale(range(0 `H')) xlabel(0(2)`H') graphregion(color(white))  ytitle("Difference of Gini") xtitle("Quarter") legend(off) name(g3, replace)

cap mkdir Figure	
graph export "Figure/F3_Newey_West.png", replace

cap mkdir Table
export excel using "Table\T3.xlsx", firstrow(variables) replace

* save dta
keep beta se_w horizon
order horizon beta se_w

rename (beta se_w) (Wealth Wealth_sd)
save "STATA_dta\06_total.dta", replace

* (2) More control
use STATA_dta\02_LP_data.dta, clear
* cut the sample
drop if time>tq(2023,4)

tsset time
gen dGini = Gini - L.Gini
gen dlnY  = log(Y) - log(L.Y)
gen dlnPi = log(Pi) - log(L.Pi)

local p = 4

* First stage
reg R MPS L(1/`p').R L(1/`p').dGini L(1/`p').dlnY L(1/`p').dlnPi L(1/`p').EBP
test MPS

* LP-IV
gen beta = .
gen se   = .
local row = 1

forvalues h = 0/`H'{
	* dependent variable
	gen Y_dif_`h' = F`h'.Gini - L.Gini
	
	* LP-IV
	qui ivregress 2sls Y_dif_`h' L(1/`p').R L(1/`p').dGini L(1/`p').dlnY L(1/`p').dlnPi L(1/`p').EBP (R = MPS L(1/`p').R L(1/`p').dGini L(1/`p').dlnY L(1/`p').dlnPi L(1/`p').EBP), r
	
	* results
	replace beta = _b[R] in `row'
	replace se   = _se[R] in `row'
	
	local ++row
}

keep beta se

gen horizon=0-1+_n
drop if horizon>`H'
gen LB = beta-se
gen UB = beta+se
gen LB2 = beta-1.645*se
gen UB2 = beta+1.645*se

tsset horizon

twoway (rarea UB LB horizon, fcolor(red%30) lcolor(%0)) (rarea UB2 LB2 horizon, fcolor(red%10) lcolor(%0)) (rarea UB2 LB2 horizon, fcolor(red%10) lcolor(%0)) (line beta horizon, lcolor(red) lwidth(thick)), yline(0, lcolor(black) lpattern(dash)) xscale(range(0 `H')) xlabel(0(2)`H') graphregion(color(white))  ytitle("Difference of Gini") xtitle("Quarter") legend(off) name(g2, replace)

cap mkdir Figure	
graph export "Figure/F4.png", replace

cap mkdir Table
export excel using "Table\T4.xlsx", firstrow(variables) replace

graph combine g1 g2, ycommon col(2) iscale(1) name(g_combined, replace) graphregion(color(white) margin(zero)) plotregion(margin(zero) lcolor(white)) xsize(10) ysize(4)

graph export "Figure/F3_F4.png", replace




