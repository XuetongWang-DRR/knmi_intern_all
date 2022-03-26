#qrf

final = read.table("E:/output/48hour_forecast_alldata.txt", sep = "\t", header = T)
final = na.omit(final)
final$R4final1 = as.numeric(final$R4final1)
library(gamlss)
library(dplyr)
final2018 = final[1:164126,]
final2020 = setdiff(final,final2018)
set.seed(1)
final2020 <- final2020 %>%
  group_by(EC4time1) %>%
  sample_n(100)

#set.seed(1)
#final2018 <- final2018 %>%
#  group_by(EC4time1) %>%
#  sample_n(50)

library(quantregForest)
library(lubridate)
library(scoringRules)
prob <- 1:51 / 52
library(verification)

#qts_save <- matrix(NA, 
#                   nrow = nrow(final2018), 
#                   ncol = length(prob))

rf <- ranger(R4final1 ~ ., final2020[,1:99], quantreg = TRUE)


qrF_model <- quantregForest(x = final2020[, !(names(final2020) %in% c("R4final1", "EC4time1", "EC4ID1","R4ID1"))], 
                            y = final2020$R4final1,
                            nodesize = 5,
                            mtry = 3,
                            ntree = 1000)
saveRDS(object = qrF_model, file = "24h_forecast_100sample_5_5_500_qrf.rds")

mod1 = readRDS(file = "24h_forecast_100sample_5_5_500_qrf.rds")
qrF_model = mod1

qrF_prediction <-   predict(qrF_model,
                            newdata = final2018[, !(names(final2018) %in% c("R4final1", "EC4time1", "EC4ID1","R4ID1"))],
                            what = prob,
                            all = TRUE)
qrF_prediction = as.data.frame(qrF_prediction)
colnames(qrF_prediction) = c(1:51)
qrF_prediction$observation = final2018$R4final1

qrf_crps <- crps_sample(y = final2018$R4final1, dat = qrF_prediction)
mean(qrf_crps)

crps3 = qrF_prediction
crps3$crps=apply(crps3,1,function(x) crps_sample(as.numeric(x['observation']),as.numeric(x[1:51]))) 
clim=crps3$observation
crps3$crps_clim=apply(crps3,1,function(x) crps_sample(as.numeric(x['observation']),clim))
crps3$crpss=(crps3$crps-crps3$crps_clim)/(-crps3$crps_clim)
crps3 = crps3[,53:55]
write.table(crps3, "crps_24h_ALLsample_5_5_500.txt", sep = "\t", col.names = T, row.names = F)

condEcdf <-   predict(qrF_model,
                            newdata = final2018[, !(names(final2018) %in% c("R4final1", "EC4time1", "EC4ID1","R4ID1"))],
                            what = ecdf,
                            all = TRUE)

CDF = as.data.frame(matrix(NA,nrow = dim(final2018)[1], ncol = 11))
for (i in 1:164126){
  CDF[i,1] = condEcdf[[i]](0.05)
  CDF[i,2] = condEcdf[[i]](0.1)
  CDF[i,3] = condEcdf[[i]](0.2)
  CDF[i,4] = condEcdf[[i]](0.5)
  CDF[i,5] = condEcdf[[i]](0.8)
  CDF[i,6] = condEcdf[[i]](1)
  CDF[i,7] = condEcdf[[i]](2)
  CDF[i,8] = condEcdf[[i]](3)
  CDF[i,9] = condEcdf[[i]](5)
  CDF[i,10] = condEcdf[[i]](7.5)
  CDF[i,11] = condEcdf[[i]](10)
}  
write.table(CDF, "CDF_48hour_100sample_5_5_500_qrf.txt", sep = "\t", col.names = T, row.names = F)
library(verification)

rm(CDF)
CDF = read.table("E:/output/CDF_48hour_100sample_5_5_500_qrf.txt", sep = "\t", header = T)
CDF = 1-CDF
#final2018 = final[1:164126,]
obs = as.data.frame(final2018$R4final1)
obs$`final2018$R4final1`[which(obs$`final2018$R4final1` <= 7.5)] <- 100
obs$`final2018$R4final1`[which(obs$`final2018$R4final1` >7.5 & obs$`final2018$R4final1` <100)] <- 1
obs$`final2018$R4final1`[which(obs$`final2018$R4final1` == 100)] <- 0
obs$`final2018$R4final1`[which(obs$`final2018$R4final1` <= 0.05)] <- 0
obs$`final2018$R4final1`[which(obs$`final2018$R4final1` > 0.05)] <- 1

df = cbind.data.frame(CDF$V1,obs$`final2018$R4final1`)
colnames(df) = c("pred","obs")
A<- verify(df$obs, df$pred, frcst.type = "prob", obs.type = "binary")

reliability.plot(A, titl = "Reliability diagram of qrf with 7.5mm threshold (48h forecast)")



  brier = brier(df$obs,df$pred, bins=F)
  brier.ss = brier$ss
  brier.bs = brier$bs



#require the functions from gbex package
setwd("D:/ITC/KNMI/KNMI/knmi_intern_all/gbex-master/R")
source('gbex.R')
