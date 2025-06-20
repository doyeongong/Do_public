* Step 04 SCF Gini
* Yeongwoong Do (2025)

clear
cls

cd "C:\Users\82106\Documents\5. 2025\MP_Wealth_inequality\replication_package"


********************************************************************************

* (1) Create the excel file
cap mkdir "Table"
putexcel set "Table/T2.xlsx", replace

* (2) Header
putexcel A1 = "year", bold
putexcel B1 = "SCF_Gini", bold
putexcel C1 = "DFA_Gini", bold

********************************************************************************

* YEAR list
local years 1989 1992 1995 1998 2001 2004 2007 2010 2013 2016 2019 2022
local row = 2

foreach y of local years {
	
	* open the file
	local short = substr("`y'", 3, 2)
    local file = "raw_data\scf\p`short'i6.dta"
	use "`file'", clear

	* Upper -> Lower
	rename _all, lower
	
	* generate YEAR
	gen YEAR = `y'
	putexcel A`row' = `y'


*************************<Weigth> **********************************************
*   sample and weight adjustments;

if YEAR== 1989 {
    gen ID = x1
	gen IID = xx1
}
else{
    gen ID = y1
	gen IID = yy1
}

*   only keep observations with valid ID and weight
keep if ID>0 & IID>0 & x42001>0

*   retain original weight: WGT0
gen WGT0=x42001
	
*************************< FIN; start >*****************************************

********************************************************************************
* 1) Liquid assets
********************************************************************************

*  1-1) checking accounts other than money market
gen CHECKING = 0
replace CHECKING = CHECKING + max(0, x3506) if x3507 == 5
replace CHECKING = CHECKING + max(0, x3510) if x3511 == 5
replace CHECKING = CHECKING + max(0, x3514) if x3515 == 5
replace CHECKING = CHECKING + max(0, x3518) if x3519 == 5
replace CHECKING = CHECKING + max(0, x3522) if x3523 == 5
replace CHECKING = CHECKING + max(0, x3526) if x3527 == 5
replace CHECKING = CHECKING + max(0, x3529) if x3527 == 5

* 1-2) prepaid cards
gen PREPAID = 0
if YEAR>=2016{
	replace PREPAID = max(0, x7596)
}

* 1-3) savings accounts
gen SAVING = 0
if YEAR <=2001{
	replace SAVING = max(0, x3804) + max(0, x3807) + max(0, x3810) + max(0, x3813) + max(0, x3816) + max(0, x3818)
}
else{
	replace SAVING = max(0, x3730) * (x3732 != 4 & x3732 != 30) + max(0, x3736) * (x3738 != 4 & x3738 != 30) + max(0, x3742) * (x3744 != 4 & x3744 != 30) + max(0, x3748) * (x3750 != 4 & x3750 != 30) + max(0, x3754) * (x3756 != 4 & x3756 != 30) + max(0, x3760) * (x3762 != 4 & x3762 != 30) + max(0, x3765) * (x3762 != 4 & x3762 != 30)
}


* 1-4)   money market deposit accounts
*   NOTE: includes money market accounts used for checking and other
*    money market account held at commercial banks, savings and
*    loans, savings banks, and credit unions;
gen MMDA = 0

if YEAR <= 2001{
	replace MMDA = max(0, x3506) * (x3507 == 1 & inrange(x9113, 11, 13)) + max(0, x3510) * (x3511 == 1 & inrange(x9114, 11, 13)) + max(0, x3514) * (x3515 == 1 & inrange(x9115, 11, 13)) + max(0, x3518) * (x3519 == 1 & inrange(x9116, 11, 13)) + max(0, x3522) * (x3523 == 1 & inrange(x9117, 11, 13)) + max(0, x3526) * (x3527 == 1 & inrange(x9118, 11, 13)) + max(0, x3529) * (x3527 == 1 & inrange(x9118, 11, 13)) + max(0, x3706) * inrange(x9131, 11, 13) + max(0, x3711) * inrange(x9132, 11, 13) + max(0, x3716) * inrange(x9133, 11, 13) + max(0, x3718) * inrange(x9133, 11, 13) 
}
else{
	replace MMDA = max(0, x3506) * (x3507 == 1 & inrange(x9113, 11, 13)) + max(0, x3510) * (x3511 == 1 & inrange(x9114, 11, 13)) + max(0, x3514) * (x3515 == 1 & inrange(x9115, 11, 13)) + max(0, x3518) * (x3519 == 1 & inrange(x9116, 11, 13)) + max(0, x3522) * (x3523 == 1 & inrange(x9117, 11, 13)) + max(0, x3526) * (x3527 == 1 & inrange(x9118, 11, 13)) + max(0, x3529) * (x3527 == 1 & inrange(x9118, 11, 13)) + max(0, x3730) * (inlist(x3732, 4, 30) & inrange(x9259, 11, 13)) + max(0, x3736) * (inlist(x3738, 4, 30) & inrange(x9260, 11, 13)) + max(0, x3742) * (inlist(x3744, 4, 30) & inrange(x9261, 11, 13)) + max(0, x3748) * (inlist(x3750, 4, 30) & inrange(x9262, 11, 13)) + max(0, x3754) * (inlist(x3756, 4, 30) & inrange(x9263, 11, 13)) + max(0, x3760) * (inlist(x3762, 4, 30) & inrange(x9264, 11, 13)) + max(0, x3765) * (inlist(x3762, 4, 30) & inrange(x9264, 11, 13)) 
}

	
* 1-5)   money market mutual funds;
*   NOTE: includes money market accounts used for checking and other
*    money market account held at institutions other than commercial
*    banks, savings and loans, savings banks, and credit unions;
gen MMMF = 0
if YEAR <=2001{
	replace MMMF = max(0, x3506)*(x3507 == 1 & (x9113 < 11 | x9113 > 13)) + max(0, x3510)*(x3511 == 1 & (x9114 < 11 | x9114 > 13)) + max(0, x3514)*(x3515 == 1 & (x9115 < 11 | x9115 > 13)) + max(0, x3518)*(x3519 == 1 & (x9116 < 11 | x9116 > 13)) + max(0, x3522)*(x3523 == 1 & (x9117 < 11 | x9117 > 13)) + max(0, x3526)*(x3527 == 1 & (x9118 < 11 | x9118 > 13)) + max(0, x3529)*(x3527 == 1 & (x9118 < 11 | x9118 > 13)) + max(0, x3706)*(x9131 < 11 | x9131 > 13) + max(0, x3711)*(x9132 < 11 | x9132 > 13) + max(0, x3716)*(x9133 < 11 | x9133 > 13) + max(0, x3718)*(x9133 < 11 | x9133 > 13) 
}
else{
	replace MMMF = max(0, x3506)*(x3507 == 1 & (x9113 < 11 | x9113 > 13)) + max(0, x3510)*(x3511 == 1 & (x9114 < 11 | x9114 > 13)) + max(0, x3514)*(x3515 == 1 & (x9115 < 11 | x9115 > 13)) + max(0, x3518)*(x3519 == 1 & (x9116 < 11 | x9116 > 13)) + max(0, x3522)*(x3523 == 1 & (x9117 < 11 | x9117 > 13)) + max(0, x3526)*(x3527 == 1 & (x9118 < 11 | x9118 > 13)) + max(0, x3529)*(x3527 == 1 & (x9118 < 11 | x9118 > 13)) + max(0, x3730)*(inlist(x3732, 4, 30) & (x9259 < 11 | x9259 > 13)) + max(0, x3736)*(inlist(x3738, 4, 30) & (x9260 < 11 | x9260 > 13)) + max(0, x3742)*(inlist(x3744, 4, 30) & (x9261 < 11 | x9261 > 13)) + max(0, x3748)*(inlist(x3750, 4, 30) & (x9262 < 11 | x9262 > 13)) + max(0, x3754)*(inlist(x3756, 4, 30) & (x9263 < 11 | x9263 > 13)) + max(0, x3760)*(inlist(x3762, 4, 30) & (x9264 < 11 | x9264 > 13)) + max(0, x3765)*(inlist(x3762, 4, 30) & (x9264 < 11 | x9264 > 13))
}
	
*  all types of money market accounts;
gen MMA = MMDA + MMMF

* 1-6) call accounts at brokerages;
gen CALL= max(0,x3930)	

* SUMMATION
gen LIQ = CHECKING+SAVING+MMA+CALL+PREPAID	

********************************************************************************	
* 2) certificates of deposit  	
********************************************************************************

* 2-a) certificates of deposit   
gen CDS=max(0,x3721)

********************************************************************************
	
********************************************************************************	
* 3) total directly-held mutual funds, excluding MMMFs
********************************************************************************

* 3-1) stock mutual funds
gen STMUTF=(x3821==1)*max(0,x3822)

* 3-2) tax-free bond mutual funds
gen TFBMUTF=(x3823==1)*max(0,x3824)

* 3-3) government bond mutual funds
gen GBMUTF=(x3825==1)*max(0,x3826)

* 3-4) other bond mutual funds
gen OBMUTF=(x3827==1)*max(0,x3828)

* 3-5) combination and other mutual funds
gen COMUTF = 0
if YEAR >= 2004{
	replace COMUTF = max(0, x3830) if x3829 == 1
}

* 3-6) other mutual funds
gen OMUTF = 0
if YEAR>=2004{
	replace OMUTF = max(0, x7787) if x7785 == 1
}


* SUMMATION
gen NMMF=STMUTF+TFBMUTF+GBMUTF+OBMUTF+COMUTF+OMUTF

********************************************************************************	
* 4) Stocks
********************************************************************************	

* stocks
gen STOCKS = max(0, x3915)

********************************************************************************		

********************************************************************************	
* 5) bonds, not including bond funds or savings bonds	
********************************************************************************

* 5-1) tax-exempt bonds (state and local bonds)
gen NOTxBND=x3910

* 5-2) mortgage-backed bonds;
gen MORTBND=x3906

* 5-3) US government and government agency bonds and bills;
gen GOVTBND=x3908

* 5-4) corporate and foreign bonds;
gen OBND = 0
if YEAR < 1992{
	replace OBND = x3912
}
	else{
	replace OBND = x7634 + x7633
}

* SUMMATION
gen BOND=NOTxBND+MORTBND+GOVTBND+OBND

********************************************************************************
* 6) total quasi-liquid: sum of IRAs, thrift accounts, and future pensions;	
********************************************************************************
* exclude RETQLIQ 
********************************************************************************


********************************************************************************
* 7) savings bonds
********************************************************************************
gen SAVBND = x3902
********************************************************************************
	
********************************************************************************
* 8) cash value of whole life insurance
********************************************************************************	
gen CASHLI=max(0,x4006)
********************************************************************************

********************************************************************************
* 9) other managed assets (trusts, annuities and managed investment accounts 
*  in which HH has equity interest)
********************************************************************************	
gen OTHMA = .

if YEAR >= 2004{
	replace OTHMA = max(0, x6577) + max(0, x6587)
}
if YEAR >= 1998 & YEAR <=2001{
	replace OTHMA = max(0, x6820)  + max(0, x6835) 
}
if YEAR<1998 | (YEAR>2001 & YEAR<2004){
	replace OTHMA = max(0, x3942) 
}
********************************************************************************

********************************************************************************
* 10) other financial assets: includes loans from the household to
*    someone else, future proceeds, royalties, futures, non-public
*    stock, deferred compensation, oil/gas/mineral invest., cash
*    n.e.c.;
********************************************************************************

* assume that PUBLIC==YES
gen OTHFIN = x4018
replace OTHFIN = OTHFIN + x4022 if inlist(x4020, 61,62,63,64,65,66,71,72,73,74,77,80,81,-7)
replace OTHFIN = OTHFIN + x4026 if inlist(x4024, 61,62,63,64,65,66,71,72,73,74,77,80,81,-7)
replace OTHFIN = OTHFIN + x4030 if inlist(x4028, 61,62,63,64,65,66,71,72,73,74,77,80,81,-7)

********************************************************************************
* SUMMATION
gen FIN = LIQ+CDS+NMMF+STOCKS+BOND+SAVBND+CASHLI+OTHMA

*************************< FIN; end >*******************************************





*************************< NFIN; start>*****************************************

* 1) value of all vehicles 
if YEAR >= 1995{
    gen VEHIC = max(0,x8166)+max(0,x8167)+max(0,x8168)+max(0,x8188)+ max(0,x2422)+max(0,x2506)+max(0,x2606)+max(0,x2623)
}
if YEAR < 1995{
    gen VEHIC=max(0,x8166)+max(0,x8167)+max(0,x8168)+max(0,x2422)+max(0,x2506)+max(0,x2606)+max(0,x2623)
}

********************************************************************************

* 2)  primary residence
* 1. x507 capped at 9000 (90%);
replace x507 = 9000 if x507 > 9000

* 2. Initialize FARMBUS = 0 and compute share ratio if applicable
gen FARMBUS = 0
gen SHARE = x507 / 10000 if x507 > 0
gen REMAIN = (10000 - x507) / 10000 if x507 > 0

* 3. Compute initial FARMBUS if x507 > 0
replace FARMBUS = SHARE * (x513 + x526 - x805 - x905 - x1005) if x507 > 0

* 4. Adjust real estate portion of farm assets
foreach var in x805 x808 x813 x905 x908 x913 x1005 x1008 x1013 {
    replace `var' = `var' * REMAIN if x507 > 0
}

* 5. Adjust based on x1103
replace FARMBUS = FARMBUS - x1108 * SHARE if x507 > 0 & x1103 == 1
replace x1108 = x1108 * REMAIN if x507 > 0 & x1103 == 1
replace x1109 = x1109 * REMAIN if x507 > 0 & x1103 == 1

* 6. Adjust based on x1114
replace FARMBUS = FARMBUS - x1119 * SHARE if x507 > 0 & x1114 == 1
replace x1119 = x1119 * REMAIN if x507 > 0 & x1114 == 1
replace x1120 = x1120 * REMAIN if x507 > 0 & x1114 == 1

* 7. Adjust based on x1125
replace FARMBUS = FARMBUS - x1130 * SHARE if x507 > 0 & x1125 == 1
replace x1130 = x1130 * REMAIN if x507 > 0 & x1125 == 1
replace x1131 = x1131 * REMAIN if x507 > 0 & x1125 == 1

* 8. Adjustment for x1136 if farm equity components exist
gen TEMP_DENOM = x1108 + x1119 + x1130
gen TEMP_NUM = x1108*(x1103==1) + x1119*(x1114==1) + x1130*(x1125==1)
gen WEIGHT = TEMP_NUM / TEMP_DENOM if TEMP_DENOM > 0

replace FARMBUS = FARMBUS - x1136 * SHARE * WEIGHT if x507 > 0 & x1136 > 0 & TEMP_DENOM > 0
replace x1136 = x1136 * REMAIN * WEIGHT if x507 > 0 & x1136 > 0 & TEMP_DENOM > 0

* 9. Compute value of primary residence
gen HOUSES = max(0,x604) + max(0,x614) + max(0,x623) + max(0,x716) + ((10000 - max(0,x507)) / 10000) * (x513 + x526)


********************************************************************************

*  3) other residential real estate: includes land contracts/notes household has made, properties other than the principal residence that are coded as 1-4 family residences, time shares, and vacations homes

gen ORESRE = .

if YEAR >=2013 {
    replace ORESRE = max(x1306, x1310) + max(x1325, x1329) + max(0, x1339) + inlist(x1703,12,14,21,22,25,40,41,42,43,44,49,50,52,999)*max(0,x1706)*(x1705/10000) + inlist(x1803,12,14,21,22,25,40,41,42,43,44,49,50,52,999)*max(0,x1806)*(x1805/10000) + max(0, x2002)
}
if YEAR == 2010{
    replace ORESRE = max(x1405, x1409) + max(x1505, x1509) + max(0, x1619) + inlist(x1703,12,14,21,22,25,40,41,42,43,44,49,50,52,999)*max(0,x1706)*(x1705/10000) + inlist(x1803,12,14,21,22,25,40,41,42,43,44,49,50,52,999)*max(0,x1806)*(x1805/10000) + max(0, x2002)
}
if (YEAR > 2010 & YEAR<2013) | (YEAR <2010){
    replace ORESRE = max(x1405, x1409) + max(x1505, x1509) + max(x1605, x1609) + max(0, x1619) + inlist(x1703,12,14,21,22,25,40,41,42,43,44,49,50,52,999)*max(0,x1706)*(x1705/10000) + inlist(x1803,12,14,21,22,25,40,41,42,43,44,49,50,52,999)*max(0,x1806)*(x1805/10000) + inlist(x1903,12,14,21,22,25,40,41,42,43,44,49,50,52,999)*max(0,x1906)*(x1905/10000) + max(0, x2002)
}


********************************************************************************

* 4) net equity in nonresidential real estate

gen NNRESRE = 0

if YEAR >= 2010 {
        replace NNRESRE = inlist(x1703, 1,2,3,4,5,6,7,10,11,13,15,24,45,46,47,48,51,53,-7) * max(0, x1706) * (x1705/10000) + inlist(x1803, 1,2,3,4,5,6,7,10,11,13,15,24,45,46,47,48,51,53,-7) * max(0, x1806) * (x1805/10000) + max(0, x2012) - inlist(x1703, 1,2,3,4,5,6,7,10,11,13,15,24,45,46,47,48,51,53,-7) * x1715 * (x1705/10000) - inlist(x1803, 1,2,3,4,5,6,7,10,11,13,15,24,45,46,47,48,51,53,-7) * x1815 * (x1805/10000) - x2016
}
    else {
        replace NNRESRE = inlist(x1703, 1,2,3,4,5,6,7,10,11,13,15,24,45,46,47,48,51,53,-7) * max(0, x1706) * (x1705/10000) + inlist(x1803, 1,2,3,4,5,6,7,10,11,13,15,24,45,46,47,48,51,53,-7) * max(0, x1806) * (x1805/10000) + inlist(x1903, 1,2,3,4,5,6,7,10,11,13,15,24,45,46,47,48,51,53,-7) * max(0, x1906) * (x1905/10000) + max(0, x2012) - inlist(x1703, 1,2,3,4,5,6,7,10,11,13,15,24,45,46,47,48,51,53,-7) * x1715 * (x1705/10000) - inlist(x1803, 1,2,3,4,5,6,7,10,11,13,15,24,45,46,47,48,51,53,-7) * x1815 * (x1805/10000) - inlist(x1903, 1,2,3,4,5,6,7,10,11,13,15,24,45,46,47,48,51,53,-7) * x1915 * (x1905/10000) - x2016
}


********************************************************************************
*  5) business interests

gen farmbus = 0

quietly {
    if x507 > 0 {
        replace farmbus = (x507/10000) * (x513 + x526 - x805 - x905 - x1005)

        local ratio = (10000 - x507)/10000

        replace x805  = x805  * `ratio'
        replace x808  = x808  * `ratio'
        replace x813  = x813  * `ratio'
        replace x905  = x905  * `ratio'
        replace x908  = x908  * `ratio'
        replace x913  = x913  * `ratio'
        replace x1005 = x1005 * `ratio'
        replace x1008 = x1008 * `ratio'
        replace x1013 = x1013 * `ratio'

        if x1103 == 1 {
            replace farmbus = farmbus - x1108 * (x507/10000)
            replace x1108 = x1108 * `ratio'
            replace x1109 = x1109 * `ratio'
        }

        if x1114 == 1 {
            replace farmbus = farmbus - x1119 * (x507/10000)
            replace x1119 = x1119 * `ratio'
            replace x1120 = x1120 * `ratio'
        }

        if x1125 == 1 {
            replace farmbus = farmbus - x1130 * (x507/10000)
            replace x1130 = x1130 * `ratio'
            replace x1131 = x1131 * `ratio'
        }

        if x1136 > 0 & (x1108 + x1119 + x1130 > 0) {
            gen double weight = (x1108*(x1103==1) + x1119*(x1114==1) + x1130*(x1125==1)) / (x1108 + x1119 + x1130)
            replace farmbus = farmbus - x1136 * (x507/10000) * weight
            replace x1136 = x1136 * `ratio' * weight
            drop weight
        }
    }
}


gen BUS = .

quietly {
    if YEAR >= 2010 {
        replace BUS = max(0,x3129) + max(0,x3124) - max(0,x3126)*(x3127==5) + max(0,x3121)*(inlist(x3122,1,6)) + max(0,x3229) + max(0,x3224) - max(0,x3226)*(x3227==5) + max(0,x3221)*(inlist(x3222,1,6)) + max(0,x3335) + farmbus + max(0,x3408) + max(0,x3412) + max(0,x3416) + max(0,x3420) + max(0,x3452) + max(0,x3428)
    }
    else {
        replace BUS = max(0,x3129) + max(0,x3124) - max(0,x3126)*(x3127==5) + max(0,x3121)*(inlist(x3122,1,6)) + max(0,x3229) + max(0,x3224) - max(0,x3226)*(x3227==5) + max(0,x3221)*(inlist(x3222,1,6)) + max(0,x3329) + max(0,x3324) - max(0,x3326)*(x3327==5) + max(0,x3321)*(inlist(x3322,1,6)) + max(0,x3335) + farmbus + max(0,x3408) + max(0,x3412) + max(0,x3416) + max(0,x3420) + max(0,x3424) + max(0,x3428)
    }
}


********************************************************************************
*  5) other nonfinancial assets

gen OTHNFIN=x4022+x4026+x4030-OTHFIN+x4018


********************************************************************************

	* Summation
	gen NFIN=VEHIC+HOUSES+ORESRE+NNRESRE+BUS+OTHNFIN
	gen DASSET=FIN+NFIN

	ineqdeco DASSET [aw=WGT0]
	local gini = r(gini)
	putexcel B`row' = `gini'

	
	* DFA_gini
	
	use STATA_dta\01_DFA_gini.dta, clear
	local timeval = tq(`y',4)
	qui summ ps_gini_DA if time == `timeval'
	local value = r(mean)	  
	putexcel C`row' = `value'
	
	local ++row

}
