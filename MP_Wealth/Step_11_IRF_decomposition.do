* Step_11_IRF_decomposition
* Yeongwoong Do (2025)

* This do-file must be run after "Step_01" and "Step_10".

clear
cls

cd "C:\Users\82106\Documents\5. 2025\MP_Wealth_inequality\replication_package"

* (1) Long-run average shares and component ginis
use STATA_dta\02_LP_data.dta, clear

keep time s_* ps_gini_*

collapse (mean) s_* ps_gini_* (sd) s_RE_sd= s_RE s_CD_sd= s_CD s_FI_sd= s_FI s_FE_sd= s_FE G_RE_sd= ps_gini_RE G_CD_sd= ps_gini_CD G_FI_sd= ps_gini_FI G_FE_sd= ps_gini_FE

gen horizon = 0
rename (ps_gini_RE ps_gini_CD ps_gini_FI ps_gini_FE) (G_RE G_CD G_FI G_FE)
order horizon

merge 1:m horizon using STATA_dta\05_share_component.dta
drop _merge

merge 1:1 horizon using STATA_dta\06_total.dta
drop _merge

* (1-b) contribution and sd by delta method
foreach var in RE CD FI FE{
	gen con_IRF_s_`var' = G_`var'[1]*beta_s_`var'
	gen con_IRF_G_`var' = s_`var'[1]*beta_G_`var'
	gen con_IRF_s_`var'_sd = sqrt(G_`var'[1]^2*se_s_`var'^2 + G_`var'_sd[1]^2*beta_s_`var'^2)
	gen con_IRF_G_`var'_sd = sqrt(s_`var'[1]^2*se_G_`var'^2 + s_`var'_sd[1]^2*beta_G_`var'^2)
}


* (2) point estimates(graph)


gen sum_in = con_IRF_s_RE + con_IRF_s_CD + con_IRF_s_FI + con_IRF_s_FE + con_IRF_G_RE + con_IRF_G_CD + con_IRF_G_FI + con_IRF_G_FE

* graph
local H = 12

label variable Wealth "total wealth IRF"
label variable con_IRF_s_FI "Interest-related FA share"
label variable con_IRF_s_FE "Equity-related FA share"

gen chart_unexp = Wealth-con_IRF_s_FE
label variable chart_unexp "others"


twoway (bar chart_unexp horizon, color(gray%50)  barwidth(0.7)) (bar con_IRF_s_FI horizon , color(blue%50)  barwidth(0.7)) (bar con_IRF_s_FE horizon, color(red%50) barwidth(0.7)) (line Wealth horizon, lcolor(cranberry) lwidth(thick)), graphregion(color(white)) yline(0, lcolor(black)) xscale(range(0 `H')) xlabel(0(2)`H') legend(size(small)) ytitle("Contribution to the total Wealth IRF") xtitle("Quarter")

cap mkdir Figure	
graph export "Figure/F7_decomposition.png", replace

* (3) Table (point)

gen share = con_IRF_s_RE + con_IRF_s_CD + con_IRF_s_FI + con_IRF_s_FE
gen compon = con_IRF_G_RE + con_IRF_G_CD + con_IRF_G_FI + con_IRF_G_FE
gen unexp  = Wealth - sum_in

cap mkdir "Table"
putexcel set "Table/T8_decomposition.xlsx", replace

* Header
putexcel B1:F1 = "Horizon(quarter)", merge
putexcel B2 = "0"
putexcel C2 = "4"
putexcel D2 = "8"
putexcel E2 = "12"
putexcel F2 = "peak(7)"

putexcel A3 = "Total wealth Gini"
putexcel A4 = "Share effects"
putexcel A5 = "    Real estate"
putexcel A6 = "    Consumer durables"
putexcel A7 = "    Financial assets(interest-related)"
putexcel A8 = "    Financial assets(equity-related)"
putexcel A9 = "Component Gini effects"
putexcel A10 = "    Real estate"
putexcel A11 = "    Consumer durables"
putexcel A12 = "    Financial assets(interest-related)"
putexcel A13 = "    Financial assets(equity-related)"
putexcel A14 = "Residuals"

* Values
* Note: horizon 0 = row 1

local row = 3
foreach varname in Wealth share con_IRF_s_RE con_IRF_s_CD con_IRF_s_FI con_IRF_s_FE compon con_IRF_G_RE con_IRF_G_CD con_IRF_G_FI con_IRF_G_FE unexp{
	putexcel B`row' = `varname'[1]
	putexcel C`row' = `varname'[5]
	putexcel D`row' = `varname'[9]
	putexcel E`row' = `varname'[13]
	putexcel F`row' = `varname'[8]
	
	local ++row
}


* (4) Table (point with SE)

* delta method


cap mkdir "Table"
putexcel set "Table/T8_decomposition_with_std.xlsx", replace

* Header
putexcel B1:F1 = "Horizon(quarter)", merge
putexcel B2 = "0"
putexcel C2 = "4"
putexcel D2 = "8"
putexcel E2 = "12"
putexcel F2 = "peak(7)"

putexcel A3 = "Total wealth Gini"
putexcel A5 = "Share effects"
putexcel A7 = "    Real estate"
putexcel A9 = "    Consumer durables"
putexcel A11 = "    Financial assets(interest-related)"
putexcel A13 = "    Financial assets(equity-related)"
putexcel A15 = "Component Gini effects"
putexcel A17 = "    Real estate"
putexcel A19 = "    Consumer durables"
putexcel A21 = "    Financial assets(interest-related)"
putexcel A23 = "    Financial assets(equity-related)"
putexcel A25 = "Residuals"

* Values
* Note: horizon 0 = row 1

local row = 3
foreach varname in Wealth share con_IRF_s_RE con_IRF_s_CD con_IRF_s_FI con_IRF_s_FE compon con_IRF_G_RE con_IRF_G_CD con_IRF_G_FI con_IRF_G_FE unexp{
	
	local strval = string(`varname'[1], "%5.4f")
	putexcel B`row' = "`strval'"
	
	local strval = string(`varname'[5], "%5.4f")	
	putexcel C`row' = "`strval'"
	
	local strval = string(`varname'[9], "%5.4f")	
	putexcel D`row' = "`strval'"
	
	local strval = string(`varname'[13], "%5.4f")	
	putexcel E`row' = "`strval'"
	
	local strval = string(`varname'[8], "%5.4f")	
	putexcel F`row' = "`strval'"
	
	local row = `row' + 2
}

local value = round(Wealth_sd[1], 0.0001)
local strval = "(" + string(`value', "%5.4f") + ")"
putexcel B4 = "`strval'"

local value = round(Wealth_sd[5], 0.0001)
local strval = "(" + string(`value', "%5.4f") + ")"
putexcel C4 = "`strval'"

local value = round(Wealth_sd[9], 0.0001)
local strval = "(" + string(`value', "%5.4f") + ")"
putexcel D4 = "`strval'"

local value = round(Wealth_sd[13], 0.0001)
local strval = "(" + string(`value', "%5.4f") + ")"
putexcel E4 = "`strval'"

local value = round(Wealth_sd[8], 0.0001)
local strval = "(" + string(`value', "%5.4f") + ")"
putexcel F4 = "`strval'"


local row = 8
foreach varname in con_IRF_s_RE_sd con_IRF_s_CD_sd con_IRF_s_FI_sd con_IRF_s_FE_sd{
	
	local strval = "(" +string(`varname'[1], "%5.4f")+ ")"
	putexcel B`row' = "`strval'"
	
	local strval = "(" +string(`varname'[5], "%5.4f")+ ")"	
	putexcel C`row' = "`strval'"
	
	local strval = "(" +string(`varname'[9], "%5.4f")+ ")"	
	putexcel D`row' = "`strval'"
	
	local strval = "(" +string(`varname'[13], "%5.4f")+ ")"	
	putexcel E`row' = "`strval'"
	
	local strval = "(" +string(`varname'[8], "%5.4f")+ ")"	
	putexcel F`row' = "`strval'"
	
	local row = `row' + 2
}


local row = 18
foreach varname in con_IRF_G_RE_sd con_IRF_G_CD_sd con_IRF_G_FI_sd con_IRF_G_FE_sd{
	
	local strval = "(" +string(`varname'[1], "%5.4f")+ ")"
	putexcel B`row' = "`strval'"
	
	local strval = "(" +string(`varname'[5], "%5.4f")+ ")"	
	putexcel C`row' = "`strval'"
	
	local strval = "(" +string(`varname'[9], "%5.4f")+ ")"	
	putexcel D`row' = "`strval'"
	
	local strval = "(" +string(`varname'[13], "%5.4f")+ ")"	
	putexcel E`row' = "`strval'"
	
	local strval = "(" +string(`varname'[8], "%5.4f")+ ")"	
	putexcel F`row' = "`strval'"
	
	local row = `row' + 2
}







