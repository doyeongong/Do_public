* Step_06_LP_IV
* Yeongwoong Do (2025)

* This do-file must be run after "Step_06".

clear
cls

cd "C:\Users\82106\Documents\5. 2025\MP_Wealth_inequality\replication_package"

****************************(Run LP-IV)*****************************************
use STATA_dta\02_LP_data.dta, clear

* cut the sample
drop if time>tq(2023,4)

tsset time

gen dlnY  = log(Y) - log(L.Y)
gen dlnPi = log(Pi) - log(L.Pi)

local p = 4 
local H = 12

rename (ps_gini_RE ps_gini_CD ps_gini_FI ps_gini_FE) (G_RE G_CD G_FI G_FE)

foreach var in RE CD FI FE{
	foreach j in s G{
	gen d`j'_`var' = `j'_`var' - L.`j'_`var'
	
		* First stage
		reg R MPS L(1/`p').R L(1/`p').d`j'_`var'
		test MPS
		
		* LP-IV
		gen beta_`j'_`var' = .
		gen se_`j'_`var' = .
		gen LB_`j'_`var' = .
		gen UB_`j'_`var' = .
		*gen LB2_`j'_`var' = .
		*gen UB2_`j'_`var' = .
		
		local row = 1
		
		forvalues h = 0/`H'{
			* dependent variable
			gen Y_dif_`j'_`var'_`h' = F`h'.`j'_`var' - L.`j'_`var'
	
			* LP-IV
			qui ivreg2 Y_dif_`j'_`var'_`h' L(1/`p').R L(1/`p').d`j'_`var' (R = MPS), r
			
			* results
			replace beta_`j'_`var' = _b[R] in `row'
			replace se_`j'_`var'   = _se[R] in `row'	
			
			local ++row
			}
		replace LB_`j'_`var'   = beta_`j'_`var'-se_`j'_`var'	
		replace UB_`j'_`var'   = beta_`j'_`var'+se_`j'_`var'
		*replace LB2_`j'_`var'   = beta_`j'_`var'-1.645*se_`j'_`var'	
		*replace UB2_`j'_`var'   = beta_`j'_`var'+1.645*se_`j'_`var'
	}
}


****************************(graph)*********************************************

keep beta* se* LB* UB*

gen horizon=0-1+_n
drop if horizon>`H'

tsset horizon

* graph by variables

foreach var in RE CD FI FE{
foreach j in s G{

	    if "`j'" == "s" {
            local ytitle "Difference of Share"
        }
        else {
            local ytitle "Difference of Gini"
        }

twoway (rarea UB_`j'_`var' LB_`j'_`var' horizon, fcolor(red%30) lcolor(%0)) (line beta_`j'_`var' horizon, lcolor(red) lwidth(thick)), yline(0, lcolor(black) lpattern(dash)) xscale(range(0 `H')) xlabel(0(2)`H') graphregion(color(white))  ytitle("`ytitle'") xtitle("Quarter") legend(off) name(g_`j'_`var', replace)

cap mkdir Figure	
graph export "Figure/F6_`j'_`var'.png", replace

}
}

* combined graph

graph combine g_s_RE g_s_CD g_s_FI g_s_FE, ycommon col(4) iscale(1) name(g_share, replace) graphregion(color(white) margin(zero)) plotregion(margin(zero) lcolor(white)) xsize(15) ysize(4)
graph export "Figure/F6_share.png", replace


graph combine g_G_RE g_G_CD g_G_FI g_G_FE, ycommon col(4) iscale(1) name(g_component_gini, replace) graphregion(color(white) margin(zero)) plotregion(margin(zero) lcolor(white)) xsize(15) ysize(4)
graph export "Figure/F6_component_Gini.png", replace


keep beta* se*

cap mkdir Table
export excel using "Table\T7.xlsx", firstrow(variables) replace

gen horizon=-1+_n
order horizon beta* se*

save "STATA_dta\05_share_component.dta", replace



