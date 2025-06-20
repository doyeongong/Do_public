* Step 01 DFA Gini
* Yeongwoong Do (2025)

clear
cls

cd "C:\Users\82106\Documents\5. 2025\MP_Wealth_inequality\replication_package"

* (1) Open the csv file
import delimited "raw_data\dfa-networth-levels-detail.csv"

* (2) Definition of variables
* RE real estate
* CD consumer durables
* FI Fixed income-related Assets (interest income)
* FE Financial Equity-related Assets
* DA Disposable Assets
rename (realestate consumerdurables) (RE CD)
gen FI = deposits + moneymarketfundshares + debtsecurities + loansassets
gen FE = corporateequitiesandmutualfundsh + miscellaneousotherequity + miscellaneousassets
gen DA = RE + CD + FI + FE

* (3) rearrange the data
keep date category RE CD FI FE DA
encode category, gen(category_num)
drop category

reshape wide RE CD FI FE DA, i(date) j(category_num)
gen time = tq(1989,2)+_n
format time %tq
drop date

foreach var in RE CD FI FE DA{
	gen total_`var' = `var'1+`var'2+`var'3+`var'4+`var'5
}

forvalues i = 1(1)5{
	foreach var in RE CD FI FE DA{
		replace `var'`i'=`var'`i'/total_`var'
	}
}

foreach var in RE CD FI FE DA{
	gen s_`var' = total_`var'/total_DA
}


drop total_*

* (4) Caculate Gini using the binned data
foreach var in RE CD FI FE DA{
	gen ps_gini_`var' = 1-2*(0.5*0.5*`var'1 + 0.4*`var'1 + 0.5*0.4*`var'2 + 0.09*(`var'1+`var'2) + 0.5*0.09*`var'3 + 0.009*(`var'1+`var'2+`var'3) + 0.5*0.009*`var'4 + 0.001*(`var'1+`var'2+`var'3+`var'4) + 0.5*0.001*`var'5)
}


keep time s_* ps_gini_*
drop s_DA

* (5) Output
cap mkdir "STATA_dta"
save "STATA_dta\01_DFA_gini.dta", replace




