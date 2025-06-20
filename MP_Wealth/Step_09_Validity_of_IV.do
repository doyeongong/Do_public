* Step_09_Validity_of_IV

* Yeongwoong Do (2025)

* This do-file must be run after "Step_01" & "Step_08".

clear
cls

cd "C:\Users\82106\Documents\5. 2025\MP_Wealth_inequality\replication_package"

**************<Generate 03_LP_data_monthly.dta>*********************************

* (1) MPS
import excel "raw_data\monetary-policy-surprises-data.xlsx", sheet("Monthly (update 2023)") firstrow

keep Year Month MPS_ORTH

save "STATA_dta\03_LP_data_monthly.dta", replace

* (2) Macrovariables

import excel "raw_data\FRED_data.xlsx", sheet("month") firstrow clear

merge m:1 Year Month using STATA_dta\03_LP_data_monthly.dta
drop _merge

* (3) RTI Gini
merge m:1 Year Month using STATA_dta\04_RTI_gini.dta
keep if _merge==3
drop _merge

gen time = tm(1988,1)-1+_n
format time %tm
order time
drop Year Month

rename (DGS1 CPIAUCSL INDPRO ebp MPS_ORTH) (R Pi Y EBP MPS)
save "STATA_dta\03_LP_data_monthly.dta", replace

cap mkdir Table
export excel using "Table\T6_dataset_month.xlsx", firstrow(variables) replace

**************<Run LP-IV>*******************************************************

use STATA_dta\03_LP_data_monthly.dta, clear


tsset time
gen dGini = Gini - L.Gini
gen dlnY  = log(Y) - log(L.Y)
gen dlnPi = log(Pi) - log(L.Pi)

local p = 12 

* First stage
reg R MPS L(1/`p').R L(1/`p').dGini
test MPS


* LP-IV
gen beta = .
gen se = .
local row = 1
local H   = 36

forvalues h = 0/`H'{
	* dependent variable
	gen Y_dif_`h' = F`h'.Gini - L.Gini
	
	* LP-IV
	qui ivreg2 Y_dif_`h' L(1/`p').R L(1/`p').dGini (R = MPS), r

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

* main
tsset horizon

twoway (rarea UB LB horizon, fcolor(red%30) lcolor(%0)) (rarea UB2 LB2 horizon, fcolor(red%10) lcolor(%0)) (line beta horizon, lcolor(red) lwidth(thick)), yline(0, lcolor(black) lpattern(dash)) xscale(range(0 `H')) xlabel(0(3)`H') graphregion(color(white))  ytitle("Difference of Gini") xtitle("Month") legend(off) name(g1, replace)

cap mkdir Figure	
graph export "Figure/F5.png", replace



