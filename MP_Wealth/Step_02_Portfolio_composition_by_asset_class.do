* Step 02 Portfolio composition by asset class
* Yeongwoong Do (2025)

clear
cls

cd "C:\Users\82106\Documents\5. 2025\MP_Wealth_inequality\replication_package"

* (1) Open the csv file
import delimited "raw_data\dfa-networth-levels-detail.csv"
keep if date=="2023:Q4"

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


* (3) share
keep date category RE CD FI FE DA

foreach var in RE CD FI FE{
	gen s_`var' = `var'/DA
}


* (4) bar chart

* category(strinng) -> cat_order(numeric)
encode category, gen(cat_order)
label define catlbl 5 "Top 0.1%" 4 "99~99.9%" 3 "90~99%" 2 "50~90%" 1 "Bottom 50%", replace
label values cat_order catlbl

graph bar s_RE s_CD s_FI s_FE, over(cat_order) stack graphregion(color(white)) ///
	ysize(3.5) xsize(5) ///
	bar(1, color(gs8%60)) ///
    bar(2, color(gs12%60)) ///
    bar(3, color(blue%60)) ///
    bar(4, color(red%60)) ///
	legend(size(3)) ///
    legend(label(1 "Real estate") ///
           label(2 "Consumer durables") ///
           label(3 "Financial asset(interest-related)") ///
           label(4 "Financial asset(equity-related)")) ///
    ylabel(0(0.2)1) ytitle("")

cap mkdir Figure	
graph export "Figure/F1.png", replace

