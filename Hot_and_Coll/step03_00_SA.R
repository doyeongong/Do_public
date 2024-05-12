### seasonal adjustment

# install.packages("forecast")
# install.packages("readxl")
# install.packages("seasonal")
# 필요한 패키지 불러오기
library(forecast)
library(readxl)
library(seasonal)

# 라이브러리 지정
setwd("C:/Users/82106/Documents/4. 2024학년도/Core_Periphery/2. Empirical/LP_HFI")

# seq() 함수를 사용하여 날짜 벡터 생성
# start_date <- as.Date("2006-01-01")
# end_date <- as.Date("2024-02-01")
# time <- seq(start_date, end_date, by = "month")

# 데이터 불러오기
ip <- read_excel("Kosis.xlsx", sheet="IP", range="B1:B291")
cpi <- read_excel("Kosis.xlsx", sheet="CPI", range="B1:B712")

# 시계열 객체로 변환
ip_ts <- ts(ip, start = c(2000, 1), frequency = 12)
cpi_ts <- ts(cpi, start = c(1965, 1), frequency = 12)
# plot.ts(cpi_ts, col='red')

# 계절성 확인
# dicom <- decompose(cpi_ts)
# dicom
# plot(dicom, col="blue")

# SA 시행
ip_sa <-seas(ip_ts)
cpi_sa <-seas(cpi_ts)

# 계열 확인
ip_sa_predict <- predict(ip_sa)
cpi_sa_predict <- predict(cpi_sa)
plot(ip_sa)


# 계절 조정된 데이터를 CSV 파일로 저장
write.csv(ip_sa, file = "ip_sa.csv", row.names = FALSE)
write.csv(cpi_sa, file = "cpi_sa.csv", row.names = FALSE)
