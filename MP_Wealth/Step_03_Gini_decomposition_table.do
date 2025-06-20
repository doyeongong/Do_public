* Step 03 Gini decomposition table
* Yeongwoong Do (2025)

* This do-file must be run after "Step_01".

clear
cls

cd "C:\Users\82106\Documents\5. 2025\MP_Wealth_inequality\replication_package"

use STATA_dta\01_DFA_gini.dta, clear

* base period = 2023Q4
local timeval = tq(2023,4)

********************************************************************************

* (1) Create the excel file
cap mkdir "Table"
putexcel set "Table/T1.xlsx", replace

* (2) Header
putexcel A1 = "Component", bold
putexcel B1 = "Share",     bold
putexcel C1 = "Within Gini", bold
putexcel D1 = "Contribution", bold

* (3) first column
putexcel A2 = "Real estate"                         
putexcel A3 = "Consumer durables"                  
putexcel A4 = "Financial asset (interest-related)" 
putexcel A5 = "Financial asset (equity-related)"    
putexcel A6 = "Total"    

* (4) other column
local row = 2
foreach var in RE CD FI FE{
    gen con_`var' = s_`var'*ps_gini_`var'
	gen con_per_`var' = con_`var'/ps_gini_DA
	
	* share
	qui summ s_`var' if time == `timeval'
    local value = round(r(mean), 0.01)
    putexcel B`row' = `value'
	
	* component Gini
	qui summ ps_gini_`var' if time == `timeval'
    local value = round(r(mean), 0.01)
    putexcel C`row' = `value'
	
	* contribution
	qui summ con_`var' if time == `timeval'
	local value = round(r(mean), 0.01)
	putexcel D`row' = `value'
	
	* contribution percent
	qui summ con_per_`var' if time == `timeval'
	local value = round(r(mean), 0.01)
	local strval = "(" + string(`value', "%3.2f") + ")"
	putexcel E`row' = "`strval'"
	
    local ++row
}

* total row
local strval = string(1.00, "%3.2f")
putexcel B6 = "`strval'"
local strval = "(" + string(1.00, "%3.2f") + ")"
putexcel E6 = "`strval'"

qui summ ps_gini_DA if time == `timeval'
local value = round(r(mean), 0.01)	  
putexcel C6 = `value'
putexcel D6 = `value'



