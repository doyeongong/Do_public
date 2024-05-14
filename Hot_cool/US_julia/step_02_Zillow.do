clear
cls

cd "C:\Users\82106\Documents\4. 2024학년도\Core_Periphery\2. Empirical\US"

* Zillow data 불러오기
import delimited Zillow_zip_code.csv, clear 

* 불필요한 변수 제거
drop sizerank regionname regiontype statename state city countyname

* metro 코드화 (인구순서로 top 30만 번호를 부여)
gen metro_id = .
replace metro_id = 1 if metro == "New York-Newark-Jersey City, NY-NJ-PA"
replace metro_id = 2 if metro == "Los Angeles-Long Beach-Anaheim, CA"
replace metro_id = 3 if metro == "Chicago-Naperville-Elgin, IL-IN-WI"
replace metro_id = 4 if metro == "Dallas-Fort Worth-Arlington, TX"
replace metro_id = 5 if metro == "Houston-The Woodlands-Sugar Land, TX"
replace metro_id = 6 if metro == "Washington-Arlington-Alexandria, DC-VA-MD-WV"
replace metro_id = 7 if metro == "Philadelphia-Camden-Wilmington, PA-NJ-DE-MD"
replace metro_id = 8 if metro == "Miami-Fort Lauderdale-Pompano Beach, FL"
replace metro_id = 9 if metro == "Atlanta-Sandy Springs-Alpharetta, GA"
replace metro_id = 10 if metro == "Boston-Cambridge-Newton, MA-NH"
replace metro_id = 11 if metro == "Phoenix-Mesa-Chandler, AZ"
replace metro_id = 12 if metro == "San Francisco-Oakland-Berkeley, CA"
replace metro_id = 13 if metro == "Riverside-San Bernardino-Ontario, CA"
replace metro_id = 14 if metro == "Detroit-Warren-Dearborn, MI"
replace metro_id = 15 if metro == "Seattle-Tacoma-Bellevue, WA"
replace metro_id = 16 if metro == "Minneapolis-St. Paul-Bloomington, MN-WI"
replace metro_id = 17 if metro == "San Diego-Chula Vista-Carlsbad, CA"
replace metro_id = 18 if metro == "Tampa-St. Petersburg-Clearwater, FL"
replace metro_id = 19 if metro == "Denver-Aurora-Lakewood, CO"
replace metro_id = 20 if metro == "Baltimore-Columbia-Towson, MD"
replace metro_id = 21 if metro == "St. Louis, MO-IL"
replace metro_id = 22 if metro == "Orlando-Kissimmee-Sanford, FL"
replace metro_id = 23 if metro == "Charlotte-Concord-Gastonia, NC-SC"
replace metro_id = 24 if metro == "San Antonio-New Braunfels, TX"
replace metro_id = 25 if metro == "Portland-Vancouver-Hillsboro, OR-WA"
replace metro_id = 26 if metro == "Pittsburgh, PA"
replace metro_id = 27 if metro == "Sacramento-Roseville-Folsom, CA"
replace metro_id = 28 if metro == "Austin-Round Rock-Georgetown, TX"
replace metro_id = 29 if metro == "Las Vegas-Henderson-Paradise, NV"
replace metro_id = 30 if metro == "Cincinnati, OH-KY-IN"

drop if metro_id==.
sort metro_id
* tab metro_id (list 조회)

* reshape
reshape long v, i(regionid metro_id) j(time)
* 결측치가 한개라도 있으면 drop
drop if v==.
by regionid, sort: gen ind=_N
keep if ind==291

* rename
rename v housing_value
* month 조정
gen month = tm(2000,1) + time - 10
format month %tm
drop time

export delimited using Zillow.csv, replace




