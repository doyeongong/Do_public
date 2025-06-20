* Step_05_DFA_vs_SCF_Gini
* Yeongwoong Do (2025)

* This do-file must be run after "Step_04".

clear
cls

cd "C:\Users\82106\Documents\5. 2025\MP_Wealth_inequality\replication_package"

import excel "Table\T2.xlsx", sheet("Sheet1") firstrow

* correlation
correlate SCF_Gini DFA_Gini

* simple linear projection
reg SCF_Gini DFA_Gini

* time series graph
tsset year

twoway ///
    (line SCF_Gini year, yaxis(1) lcolor(blue) lpattern(solid)) ///
    (line DFA_Gini year, yaxis(2) lcolor(red) lpattern(dash)), ///
    legend(label(1 "SCF Gini (left axis)") label(2 "DFA Gini (right axis)")) ///
    ytitle("SCF Gini", axis(1)) ///
    ytitle("DFA Gini", axis(2)) ///
	ylabel(0.72(0.02)0.82, axis(1)) ///
    ylabel(0.62(0.02)0.72, axis(2)) ///
    xtitle("Year") ///
    graphregion(color(white))

	
cap mkdir Figure	
graph export "Figure/F2.png", replace