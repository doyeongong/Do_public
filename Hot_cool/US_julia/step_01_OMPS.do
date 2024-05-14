* example analysis (San Diego)

clear
cls
cd "C:\Users\82106\Documents\4. 2024학년도\Core_Periphery\2. Empirical\US"


******** MP shock **************
* empty dta file
set obs 383
gen month = ym(1988, 2) + _n - 1
format month %tm
save OMPS.dta, replace


* Load OMPS(orthogonalized monetary policy surprise suggest by Bauer and Swanson(2023))
clear
import excel "FOMC_Bauer_Swanson.xlsx", sheet("High-Freq FOMC Announcemt Surp") cellrange(A1:W324) firstrow

* drop 2001-09-17
drop if Date==td(17sep2001)
keep Date MPS_ORTH
destring MPS_ORTH, replace

* duplicated list Date (it suggests there are duplicates in time variables)
collapse (sum) MPS_ORTH, by(Date)
tsset Date
* change frequency to month
gen month = ym(year(Date), month(Date))
format month %tm
collapse (sum) MPS_ORTH, by(month)

* complete dta file
merge 1:m month using OMPS.dta
replace MPS_ORTH = 0 if _merge == 2
drop _merge
sort month

export delimited using "C:\Users\82106\Documents\4. 2024학년도\교수님 RA\IPUMS\OMPS.csv", replace
* save OMPS.dta, replace 




