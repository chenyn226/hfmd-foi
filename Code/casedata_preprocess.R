
library(readxl)
library(data.table)

data <- read_xlsx("case_data_hfmd.xlsx")
data <- as.data.table(data)

ind1 <- grep("CV-A6|CV-A10|CV-A16|EV-A71",data[,`lab results`])
data1 <- data[ind1,]

ind_CVA6 <- grep("CV-A6",data1[,`lab results`])
data_CVA6 <- data1[ind_CVA6,]
data_CVA6[,virus := "CV-A6"]

ind_CVA10 <- grep("CV-A10",data1[,`lab results`])
data_CVA10 <- data1[ind_CVA10,]
data_CVA10[,virus := "CV-A10"]

ind_CVA16 <- grep("CV-A16",data1[,`lab results`])
data_CVA16 <- data1[ind_CVA16,]
data_CVA16[,virus := "CV-A16"]

ind_EVA71 <- grep("EV-A71",data1[,`lab results`])
data_EVA71 <- data1[ind_EVA71,]
data_EVA71[,virus := "EV-A71"]

data2 <- rbind(data_CVA6,data_CVA10,data_CVA16,data_EVA71)

data2[,year_sampling := as.integer(format(`sampling date`,format="%Y"))]
data2[,year_illness := as.integer(format(DateofillnessOnset,format="%Y"))]
data2[,year_end := 0]

for (i in 1:dim(data2)[1]){
  if (is.na(data2[i,DateofillnessOnset])) 
    set(data2,i,"year_end",data2[i,year_sampling])
  else if ((!is.na(data2[i,DateofillnessOnset]))&(data2[i,year_sampling <= year_illness]))
    set(data2,i,"year_end",data2[i,year_sampling])
  else
    set(data2,i,"year_end",data2[i,year_illness])
}
stopifnot(data2[,year_sampling>=year_end])

# assume the age in days/years means that at the time of sampling, not illness
# remove those with age<0, after removing these samples, no missing data
# in this process, age_integer = year_end-year_start

index_remove <- which(data2[,(is.na(`Age in years`))|(`Age in years`<0)])
data3 <- data2[-index_remove,]
data3[,birth_date := `sampling date`-lubridate::days(`Age in days`)]
data3[,year_start := as.integer(format(birth_date,format="%Y"))]
stopifnot(data3[,year_start<=year_end])
data3[,age_integer := year_end-year_start]

saveRDS(data3,"case_data_preprocessed.rds")
