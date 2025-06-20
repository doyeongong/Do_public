* Step 09 RTI(Realtimeinequality) Gini
* Yeongwoong Do (2025)

clear
cls

cd "C:\Users\82106\Documents\5. 2025\MP_Wealth_inequality\replication_package"

* (1) Open the csv file
import excel "raw_data\Wealth_Realtimeinequality.xlsx", sheet("wealth") firstrow

* (2) percentile

gen m90_99 = Top10-Top1
gen m99_99_9 = Top1-Top01
gen m99_99_99 = Top01-Top001


* make sure that summation is equal to 1
gen sum = Bottom50+ Middle40+ m90_99+ m99_99_9+ m99_99_99+ Top001
gen p00_50 = Bottom50/sum
gen p50_90 = Middle40/sum
gen p90_99 = m90_99/sum
gen p99_999 = m99_99_9/sum
gen p999_9999 = m99_99_99/sum
gen p9999_1 = Top001/sum

* (3) Caculate Gini using the binned data
gen Gini = 1-2*(0.5*0.5*p00_50 + 0.4*p00_50 + 0.5*0.4*p50_90 + 0.09*(p00_50+p50_90) + 0.5*0.09*p90_99 + 0.009*(p00_50+p50_90+p90_99) + 0.5*0.009*p99_999 + 0.001*(p00_50+p50_90+p90_99+p99_999) + 0.5*0.0009*p999_9999+0.0001*(p00_50+p50_90+p90_99+p99_999+p999_9999)+0.5*0.0001*p9999_1)

keep Year Month Gini

* (4) Output
cap mkdir "STATA_dta"
save "STATA_dta\04_RTI_gini.dta", replace





