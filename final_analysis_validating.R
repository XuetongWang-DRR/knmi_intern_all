#setwd("E:/output/final_analysis/24hours")
crps_24h_raw = read.table("E:/output/validating_model/24hours/crps_24h_raw_forecast.txt", sep = "\t", header = T)
crps_24h_zaga = read.table("E:/output/validating_model/24hours/crps_24h_3steps.txt", sep = "\t", header = T)
crps_24h_qrf = read.table("E:/output/validating_model/24hours/crps_24h_ALLsample_5_3_500.txt", sep = "\t", header = T)
CDF_24h_zaga = read.table("E:/output/validating_model/24hours/NCDF_24hour_3steps.txt", sep = "\t", header = T)
CDF_24h_qrf = read.table("E:/output/validating_model/24hours/NCDF_24hour_ALLsample_5_3_500_qrf.txt", sep = "\t", header = T)


#setwd("E:/output/final_analysis/48hours")
crps_48h_raw = read.table("E:/output/validating_model/48hours/crps_48h_raw_forecast.txt", sep = "\t", header = T)
crps_48h_zaga = read.table("E:/output/validating_model/48hours/crps_48h_3steps.txt", sep = "\t", header = T)
crps_48h_qrf = read.table("E:/output/validating_model/48hours/crps_48h_ALLsample_5_3_500.txt", sep = "\t", header = T)
CDF_48h_zaga = read.table("E:/output/validating_model/48hours/NCDF_48hour_3steps.txt", sep = "\t", header = T)
CDF_48h_qrf = read.table("E:/output/validating_model/48hours/NCDF_48hour_ALLsample_5_3_500_qrf.txt", sep = "\t", header = T)


#setwd("E:/output/final_analysis/72hours")
crps_72h_raw = read.table("E:/output/validating_model/72hours/crps_72h_raw_forecast.txt", sep = "\t", header = T)
crps_72h_zaga = read.table("E:/output/validating_model/72hours/crps_72h_3steps.txt", sep = "\t", header = T)
crps_72h_qrf = read.table("E:/output/validating_model/72hours/crps_72h_ALLsample_5_3_500.txt", sep = "\t", header = T)
CDF_72h_zaga = read.table("E:/output/validating_model/72hours/NCDF_72hour_3steps.txt", sep = "\t", header = T)
CDF_72h_qrf = read.table("E:/output/validating_model/72hours/NCDF_72hour_ALLsample_5_3_500_qrf.txt", sep = "\t", header = T)

#setwd("E:/output/final_analysis/96hours")
crps_96h_raw = read.table("E:/output/validating_model/96hours/crps_96h_raw_forecast.txt", sep = "\t", header = T)
crps_96h_zaga = read.table("E:/output/validating_model/96hours/crps_96h_3steps.txt", sep = "\t", header = T)
crps_96h_qrf = read.table("E:/output/validating_model/96hours/crps_96h_ALLsample_5_3_500.txt", sep = "\t", header = T)
CDF_96h_zaga = read.table("E:/output/validating_model/96hours/NCDF_96hour_3steps.txt", sep = "\t", header = T)
CDF_96h_qrf = read.table("E:/output/validating_model/96hours/NCDF_96hour_ALLsample_5_3_500_qrf.txt", sep = "\t", header = T)


#setwd("E:/output/final_analysis/120hours")
crps_120h_raw = read.table("E:/output/validating_model/120hours/crps_120h_raw_forecast.txt", sep = "\t", header = T)
crps_120h_zaga = read.table("E:/output/validating_model/120hours/crps_120h_3steps.txt", sep = "\t", header = T)
crps_120h_qrf = read.table("E:/output/validating_model/120hours/crps_120h_ALLsample_5_3_500.txt", sep = "\t", header = T)
CDF_120h_zaga = read.table("E:/output/validating_model/120hours/NCDF_120hour_3steps.txt", sep = "\t", header = T)
CDF_120h_qrf = read.table("E:/output/validating_model/120hours/NCDF_120hour_ALLsample_5_3_500_qrf.txt", sep = "\t", header = T)

#setwd("E:/output/final_analysis/144hours")
crps_144h_raw = read.table("E:/output/validating_model/144hours/crps_144h_raw_forecast.txt", sep = "\t", header = T)
crps_144h_zaga = read.table("E:/output/validating_model/144hours/crps_144h_3steps.txt", sep = "\t", header = T)
crps_144h_qrf = read.table("E:/output/validating_model/144hours/crps_144h_ALLsample_5_3_500.txt", sep = "\t", header = T)
CDF_144h_zaga = read.table("E:/output/validating_model/144hours/NCDF_144hour_3steps.txt", sep = "\t", header = T)
CDF_144h_qrf = read.table("E:/output/validating_model/144hours/NCDF_144hour_ALLsample_5_3_500_qrf.txt", sep = "\t", header = T)



x <- rep(c(24,48,72,96,120,144), each = 3)
y <- rep(c('raw','zaga','qrf'),times = 6)
value <- c(mean(crps_24h_raw$crpss),mean(crps_24h_zaga$crpss),mean(crps_24h_qrf$crpss),
           mean(crps_48h_raw$crpss),mean(crps_48h_zaga$crpss),mean(crps_48h_qrf$crpss),
           mean(crps_72h_raw$crpss),mean(crps_72h_zaga$crpss),mean(crps_72h_qrf$crpss),
           mean(crps_96h_raw$crpss),mean(crps_96h_zaga$crpss),mean(crps_96h_qrf$crpss),
           mean(crps_120h_raw$crpss),mean(crps_120h_zaga$crpss),mean(crps_120h_qrf$crpss),
           mean(crps_144h_raw$crpss),mean(crps_144h_zaga$crpss),mean(crps_144h_qrf$crpss))

#z <- c()
df <- data.frame(x = x, y=y, value=value)
library(ggplot2)
#不作任何条形宽度和条形距离的调整
ggplot(data = df, mapping = aes(x = factor(x), y = value, fill = y)) + geom_bar(stat = 'identity', position = 'dodge')+xlab("lead time (h)")+ylab("crpss")+
  ggtitle("crpss comparasion") +
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))

library(verification)

#24h_zaga_0.05
CDF_24h_zaga = 1- CDF_24h_zaga
final_24h = read.table("E:/output/validating_model/24hours/24hour_forecast_alldata.txt", sep = "\t", header = T)

final_24h = na.omit(final_24h)
final_24h$R4final1 = as.numeric(final_24h$R4final1)
final2018_24h = final_24h[1:164126,]
final2020_24h = setdiff(final_24h,final2018_24h)

obs = as.data.frame(final2018_24h$R4final1)
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` <= 0.05)] <- 0
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` > 0.05)] <- 1

df_24h_zaga = cbind.data.frame(CDF_24h_zaga$V1,obs$`final2018_24h$R4final1`)

colnames(df_24h_zaga) = c("pred","obs")

#A<- verify(df_24h_zaga$obs, df_24h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 0.05mm threshold (24h forecast)")


brier_24h_zaga = brier(df_24h_zaga$obs,df_24h_zaga$pred, bins=F)
brier.ss_24h_zaga_0.05 = brier_24h_zaga$ss
brier.bs_24h_zaga_0.05 = brier_24h_zaga$bs


#crps3$crps=apply(crps3,1,function(x) crps_sample(as.numeric(x['observation']),as.numeric(x[1:51]))) 
#clim=crps3$observation
#crps3$crps_clim=apply(crps3,1,function(x) crps_sample(as.numeric(x['observation']),clim))
#crps3$crpss=(mean(crps3$crps)-mean(crps3$crps_clim))/(-mean(crps3$crps_clim))
#crps3 = crps3[,53:55]


#24h_qrf_0.05
CDF_24h_qrf = 1- CDF_24h_qrf
obs = as.data.frame(final2018_24h$R4final1)
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` <= 0.05)] <- 0
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` > 0.05)] <- 1

df_24h_qrf = cbind.data.frame(CDF_24h_qrf$V1,obs$`final2018_24h$R4final1`)

colnames(df_24h_qrf) = c("pred","obs")

brier_24h_qrf = brier(df_24h_qrf$obs,df_24h_qrf$pred, bins=F)
brier.ss_24h_qrf_0.05 = brier_24h_qrf$ss
brier.bs_24h_qrf_0.05 = brier_24h_qrf$bs
rm(brier_24h_zaga,brier_24h_qrf,df_24h_zaga,df_24h_qrf,obs)
#24h_zaga_0.5

obs = as.data.frame(final2018_24h$R4final1)
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` <= 0.5)] <- 0
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` > 0.5)] <- 1

df_24h_zaga = cbind.data.frame(CDF_24h_zaga$V4,obs$`final2018_24h$R4final1`)

colnames(df_24h_zaga) = c("pred","obs")

brier_24h_zaga = brier(df_24h_zaga$obs,df_24h_zaga$pred, bins=F)
brier.ss_24h_zaga_0.5 = brier_24h_zaga$ss
brier.bs_24h_zaga_0.5 = brier_24h_zaga$bs


#24h_qrf_0.5
obs = as.data.frame(final2018_24h$R4final1)
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` <= 0.5)] <- 0
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` > 0.5)] <- 1

df_24h_qrf = cbind.data.frame(CDF_24h_qrf$V4,obs$`final2018_24h$R4final1`)

colnames(df_24h_qrf) = c("pred","obs")

brier_24h_qrf = brier(df_24h_qrf$obs,df_24h_qrf$pred, bins=F)
brier.ss_24h_qrf_0.5 = brier_24h_qrf$ss
brier.bs_24h_qrf_0.5 = brier_24h_qrf$bs

rm(brier_24h_zaga,brier_24h_qrf,df_24h_zaga,df_24h_qrf,obs)

#24h_zaga_1

obs = as.data.frame(final2018_24h$R4final1)
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` <= 1)] <- 0
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` > 1)] <- 1

df_24h_zaga = cbind.data.frame(CDF_24h_zaga$V6,obs$`final2018_24h$R4final1`)

colnames(df_24h_zaga) = c("pred","obs")

brier_24h_zaga = brier(df_24h_zaga$obs,df_24h_zaga$pred, bins=F)
brier.ss_24h_zaga_1 = brier_24h_zaga$ss

#24h_qrf_1
obs = as.data.frame(final2018_24h$R4final1)
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` <= 1)] <- 0
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` > 1)] <- 1

df_24h_qrf = cbind.data.frame(CDF_24h_qrf$V6,obs$`final2018_24h$R4final1`)

colnames(df_24h_qrf) = c("pred","obs")

brier_24h_qrf = brier(df_24h_qrf$obs,df_24h_qrf$pred, bins=F)
brier.ss_24h_qrf_1 = brier_24h_qrf$ss
rm(brier_24h_zaga,brier_24h_qrf,df_24h_zaga,df_24h_qrf,obs)


#24h_zaga_5

obs = as.data.frame(final2018_24h$R4final1)
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` <= 5)] <- 0
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` > 5)] <- 1

df_24h_zaga = cbind.data.frame(CDF_24h_zaga$V9,obs$`final2018_24h$R4final1`)

colnames(df_24h_zaga) = c("pred","obs")

#A<- verify(df_24h_zaga$obs, df_24h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 5mm threshold (24h forecast)")

brier_24h_zaga = brier(df_24h_zaga$obs,df_24h_zaga$pred, bins=F)
brier.ss_24h_zaga_5 = brier_24h_zaga$ss

#24h_qrf_5
obs = as.data.frame(final2018_24h$R4final1)
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` <= 5)] <- 0
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` > 5)] <- 1

df_24h_qrf = cbind.data.frame(CDF_24h_qrf$V9,obs$`final2018_24h$R4final1`)

colnames(df_24h_qrf) = c("pred","obs")

#A<- verify(df_24h_qrf$obs, df_24h_qrf$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of qrf with 5mm threshold (24h forecast)")


brier_24h_qrf = brier(df_24h_qrf$obs,df_24h_qrf$pred, bins=F)
brier.ss_24h_qrf_5 = brier_24h_qrf$ss
rm(brier_24h_zaga,brier_24h_qrf,df_24h_zaga,df_24h_qrf,obs)

#24h_zaga_10

obs = as.data.frame(final2018_24h$R4final1)
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` <= 10)] <- 0
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` >10)] <- 1

df_24h_zaga = cbind.data.frame(CDF_24h_zaga$V11,obs$`final2018_24h$R4final1`)

colnames(df_24h_zaga) = c("pred","obs")


#A<- verify(df_24h_zaga$obs, df_24h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 10mm threshold (24h forecast)")

brier_24h_zaga = brier(df_24h_zaga$obs,df_24h_zaga$pred, bins=F)
brier.ss_24h_zaga_10 = brier_24h_zaga$ss

#24h_qrf_10
obs = as.data.frame(final2018_24h$R4final1)
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` <=10)] <- 0
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` > 10)] <- 1

df_24h_qrf = cbind.data.frame(CDF_24h_qrf$V11,obs$`final2018_24h$R4final1`)

colnames(df_24h_qrf) = c("pred","obs")


#A<- verify(df_24h_qrf$obs, df_24h_qrf$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of qrf with 10mm threshold (24h forecast)")

brier_24h_qrf = brier(df_24h_qrf$obs,df_24h_qrf$pred, bins=F)
brier.ss_24h_qrf_10 = brier_24h_qrf$ss
rm(brier_24h_zaga,brier_24h_qrf,df_24h_zaga,df_24h_qrf,obs)



#24h_zaga_15

obs = as.data.frame(final2018_24h$R4final1)
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` <= 15)] <- 0
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` >15)] <- 1

df_24h_zaga = cbind.data.frame(CDF_24h_zaga$V13,obs$`final2018_24h$R4final1`)

colnames(df_24h_zaga) = c("pred","obs")


#A<- verify(df_24h_zaga$obs, df_24h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 15mm threshold (24h forecast)")


brier_24h_zaga = brier(df_24h_zaga$obs,df_24h_zaga$pred, bins=F)
brier.ss_24h_zaga_15 = brier_24h_zaga$ss

#24h_qrf_15
obs = as.data.frame(final2018_24h$R4final1)
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` <=15)] <- 0
obs$`final2018_24h$R4final1`[which(obs$`final2018_24h$R4final1` > 15)] <- 1

df_24h_qrf = cbind.data.frame(CDF_24h_qrf$V13,obs$`final2018_24h$R4final1`)

colnames(df_24h_qrf) = c("pred","obs")

#A<- verify(df_24h_qrf$obs, df_24h_qrf$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of qrf with 15mm threshold (24h forecast)")


brier_24h_qrf = brier(df_24h_qrf$obs,df_24h_qrf$pred, bins=F)
brier.ss_24h_qrf_15 = brier_24h_qrf$ss
rm(brier_24h_zaga,brier_24h_qrf,df_24h_zaga,df_24h_qrf,obs)

rm(final_24h,final2018_24h,final2020_24h)



#48h_zaga_0.05
CDF_48h_zaga = 1- CDF_48h_zaga
final_48h = read.table("E:/output/validating_model/48hours/48hour_forecast_alldata.txt", sep = "\t", header = T)

final_48h = na.omit(final_48h)
final_48h$R4final1 = as.numeric(final_48h$R4final1)
final2018_48h = final_48h[1:164126,]
final2020_48h = setdiff(final_48h,final2018_48h)

obs = as.data.frame(final2018_48h$R4final1)
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` <= 0.05)] <- 0
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` > 0.05)] <- 1

df_48h_zaga = cbind.data.frame(CDF_48h_zaga$V1,obs$`final2018_48h$R4final1`)

colnames(df_48h_zaga) = c("pred","obs")


#A<- verify(df_48h_zaga$obs, df_48h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 0.05mm threshold (48h forecast)")

brier_48h_zaga = brier(df_48h_zaga$obs,df_48h_zaga$pred, bins=F)
brier.ss_48h_zaga_0.05 = brier_48h_zaga$ss

#48h_qrf_0.05
CDF_48h_qrf = 1- CDF_48h_qrf
obs = as.data.frame(final2018_48h$R4final1)
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` <= 0.05)] <- 0
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` > 0.05)] <- 1

df_48h_qrf = cbind.data.frame(CDF_48h_qrf$V1,obs$`final2018_48h$R4final1`)

colnames(df_48h_qrf) = c("pred","obs")

#A<- verify(df_48h_qrf$obs, df_48h_qrf$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of qrf with 0.05mm threshold (48h forecast)")

brier_48h_qrf = brier(df_48h_qrf$obs,df_48h_qrf$pred, bins=F)
brier.ss_48h_qrf_0.05 = brier_48h_qrf$ss
rm(brier_48h_zaga,brier_48h_qrf,df_48h_zaga,df_48h_qrf,obs)

#48h_zaga_0.5

obs = as.data.frame(final2018_48h$R4final1)
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` <= 0.5)] <- 0
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` > 0.5)] <- 1

df_48h_zaga = cbind.data.frame(CDF_48h_zaga$V4,obs$`final2018_48h$R4final1`)

colnames(df_48h_zaga) = c("pred","obs")

brier_48h_zaga = brier(df_48h_zaga$obs,df_48h_zaga$pred, bins=F)
brier.ss_48h_zaga_0.5 = brier_48h_zaga$ss

#48h_qrf_0.5
obs = as.data.frame(final2018_48h$R4final1)
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` <= 0.5)] <- 0
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` > 0.5)] <- 1

df_48h_qrf = cbind.data.frame(CDF_48h_qrf$V4,obs$`final2018_48h$R4final1`)

colnames(df_48h_qrf) = c("pred","obs")

brier_48h_qrf = brier(df_48h_qrf$obs,df_48h_qrf$pred, bins=F)
brier.ss_48h_qrf_0.5 = brier_48h_qrf$ss
rm(brier_48h_zaga,brier_48h_qrf,df_48h_zaga,df_48h_qrf,obs)

#48h_zaga_1

obs = as.data.frame(final2018_48h$R4final1)
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` <= 1)] <- 0
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` > 1)] <- 1

df_48h_zaga = cbind.data.frame(CDF_48h_zaga$V6,obs$`final2018_48h$R4final1`)

colnames(df_48h_zaga) = c("pred","obs")

brier_48h_zaga = brier(df_48h_zaga$obs,df_48h_zaga$pred, bins=F)
brier.ss_48h_zaga_1 = brier_48h_zaga$ss

#48h_qrf_1
obs = as.data.frame(final2018_48h$R4final1)
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` <= 1)] <- 0
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` > 1)] <- 1

df_48h_qrf = cbind.data.frame(CDF_48h_qrf$V6,obs$`final2018_48h$R4final1`)

colnames(df_48h_qrf) = c("pred","obs")

brier_48h_qrf = brier(df_48h_qrf$obs,df_48h_qrf$pred, bins=F)
brier.ss_48h_qrf_1 = brier_48h_qrf$ss
rm(brier_48h_zaga,brier_48h_qrf,df_48h_zaga,df_48h_qrf,obs)


#48h_zaga_5

obs = as.data.frame(final2018_48h$R4final1)
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` <= 5)] <- 0
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` > 5)] <- 1

df_48h_zaga = cbind.data.frame(CDF_48h_zaga$V9,obs$`final2018_48h$R4final1`)

colnames(df_48h_zaga) = c("pred","obs")

#A<- verify(df_48h_zaga$obs, df_48h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 5mm threshold (48h forecast)")

brier_48h_zaga = brier(df_48h_zaga$obs,df_48h_zaga$pred, bins=F)
brier.ss_48h_zaga_5 = brier_48h_zaga$ss

#48h_qrf_5
obs = as.data.frame(final2018_48h$R4final1)
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` <= 5)] <- 0
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` > 5)] <- 1

df_48h_qrf = cbind.data.frame(CDF_48h_qrf$V9,obs$`final2018_48h$R4final1`)

colnames(df_48h_qrf) = c("pred","obs")
#A<- verify(df_48h_qrf$obs, df_48h_qrf$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of qrf with 5mm threshold (48h forecast)")


brier_48h_qrf = brier(df_48h_qrf$obs,df_48h_qrf$pred, bins=F)
brier.ss_48h_qrf_5 = brier_48h_qrf$ss
rm(brier_48h_zaga,brier_48h_qrf,df_48h_zaga,df_48h_qrf,obs)

#48h_zaga_10

obs = as.data.frame(final2018_48h$R4final1)
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` <= 10)] <- 0
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` >10)] <- 1

df_48h_zaga = cbind.data.frame(CDF_48h_zaga$V11,obs$`final2018_48h$R4final1`)

colnames(df_48h_zaga) = c("pred","obs")


#A<- verify(df_48h_zaga$obs, df_48h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 10mm threshold (48h forecast)")


brier_48h_zaga = brier(df_48h_zaga$obs,df_48h_zaga$pred, bins=F)
brier.ss_48h_zaga_10 = brier_48h_zaga$ss

#48h_qrf_10
obs = as.data.frame(final2018_48h$R4final1)
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` <=10)] <- 0
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` > 10)] <- 1

df_48h_qrf = cbind.data.frame(CDF_48h_qrf$V11,obs$`final2018_48h$R4final1`)

colnames(df_48h_qrf) = c("pred","obs")

#A<- verify(df_48h_qrf$obs, df_48h_qrf$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of qrf with 10mm threshold (48h forecast)")


brier_48h_qrf = brier(df_48h_qrf$obs,df_48h_qrf$pred, bins=F)
brier.ss_48h_qrf_10 = brier_48h_qrf$ss
rm(brier_48h_zaga,brier_48h_qrf,df_48h_zaga,df_48h_qrf,obs)



#48h_zaga_15

obs = as.data.frame(final2018_48h$R4final1)
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` <= 15)] <- 0
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` >15)] <- 1

df_48h_zaga = cbind.data.frame(CDF_48h_zaga$V13,obs$`final2018_48h$R4final1`)

colnames(df_48h_zaga) = c("pred","obs")


#A<- verify(df_48h_zaga$obs, df_48h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 15mm threshold (48h forecast)")

brier_48h_zaga = brier(df_48h_zaga$obs,df_48h_zaga$pred, bins=F)
brier.ss_48h_zaga_15 = brier_48h_zaga$ss

#48h_qrf_15
obs = as.data.frame(final2018_48h$R4final1)
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` <=15)] <- 0
obs$`final2018_48h$R4final1`[which(obs$`final2018_48h$R4final1` > 15)] <- 1

df_48h_qrf = cbind.data.frame(CDF_48h_qrf$V13,obs$`final2018_48h$R4final1`)

colnames(df_48h_qrf) = c("pred","obs")

#A<- verify(df_48h_qrf$obs, df_48h_qrf$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of qrf with 15mm threshold (48h forecast)")

brier_48h_qrf = brier(df_48h_qrf$obs,df_48h_qrf$pred, bins=F)
brier.ss_48h_qrf_15 = brier_48h_qrf$ss
rm(brier_48h_zaga,brier_48h_qrf,df_48h_zaga,df_48h_qrf,obs)

rm(final_48h,final2018_48h,final2020_48h)



#72h_zaga_0.05
CDF_72h_zaga = 1- CDF_72h_zaga
final_72h = read.table("E:/output/validating_model/72hours/72hour_forecast_alldata.txt", sep = "\t", header = T)

final_72h = na.omit(final_72h)
final_72h$R4final1 = as.numeric(final_72h$R4final1)
final2018_72h = final_72h[1:164126,]
final2020_72h = setdiff(final_72h,final2018_72h)

obs = as.data.frame(final2018_72h$R4final1)
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` <= 0.05)] <- 0
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` > 0.05)] <- 1

df_72h_zaga = cbind.data.frame(CDF_72h_zaga$V1,obs$`final2018_72h$R4final1`)

colnames(df_72h_zaga) = c("pred","obs")


#A<- verify(df_72h_zaga$obs, df_72h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 0.05mm threshold (72h forecast)")

brier_72h_zaga = brier(df_72h_zaga$obs,df_72h_zaga$pred, bins=F)
brier.ss_72h_zaga_0.05 = brier_72h_zaga$ss

#72h_qrf_0.05
CDF_72h_qrf = 1- CDF_72h_qrf
obs = as.data.frame(final2018_72h$R4final1)
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` <= 0.05)] <- 0
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` > 0.05)] <- 1

df_72h_qrf = cbind.data.frame(CDF_72h_qrf$V1,obs$`final2018_72h$R4final1`)

colnames(df_72h_qrf) = c("pred","obs")

#A<- verify(df_72h_qrf$obs, df_72h_qrf$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of qrf with 0.05mm threshold (72h forecast)")

brier_72h_qrf = brier(df_72h_qrf$obs,df_72h_qrf$pred, bins=F)
brier.ss_72h_qrf_0.05 = brier_72h_qrf$ss
rm(brier_72h_zaga,brier_72h_qrf,df_72h_zaga,df_72h_qrf,obs)

#72h_zaga_0.5

obs = as.data.frame(final2018_72h$R4final1)
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` <= 0.5)] <- 0
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` > 0.5)] <- 1

df_72h_zaga = cbind.data.frame(CDF_72h_zaga$V4,obs$`final2018_72h$R4final1`)

colnames(df_72h_zaga) = c("pred","obs")

brier_72h_zaga = brier(df_72h_zaga$obs,df_72h_zaga$pred, bins=F)
brier.ss_72h_zaga_0.5 = brier_72h_zaga$ss

#72h_qrf_0.5
obs = as.data.frame(final2018_72h$R4final1)
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` <= 0.5)] <- 0
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` > 0.5)] <- 1

df_72h_qrf = cbind.data.frame(CDF_72h_qrf$V4,obs$`final2018_72h$R4final1`)

colnames(df_72h_qrf) = c("pred","obs")

brier_72h_qrf = brier(df_72h_qrf$obs,df_72h_qrf$pred, bins=F)
brier.ss_72h_qrf_0.5 = brier_72h_qrf$ss
rm(brier_72h_zaga,brier_72h_qrf,df_72h_zaga,df_72h_qrf,obs)

#72h_zaga_1

obs = as.data.frame(final2018_72h$R4final1)
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` <= 1)] <- 0
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` > 1)] <- 1

df_72h_zaga = cbind.data.frame(CDF_72h_zaga$V6,obs$`final2018_72h$R4final1`)

colnames(df_72h_zaga) = c("pred","obs")

brier_72h_zaga = brier(df_72h_zaga$obs,df_72h_zaga$pred, bins=F)
brier.ss_72h_zaga_1 = brier_72h_zaga$ss

#72h_qrf_1
obs = as.data.frame(final2018_72h$R4final1)
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` <= 1)] <- 0
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` > 1)] <- 1

df_72h_qrf = cbind.data.frame(CDF_72h_qrf$V6,obs$`final2018_72h$R4final1`)

colnames(df_72h_qrf) = c("pred","obs")

brier_72h_qrf = brier(df_72h_qrf$obs,df_72h_qrf$pred, bins=F)
brier.ss_72h_qrf_1 = brier_72h_qrf$ss
rm(brier_72h_zaga,brier_72h_qrf,df_72h_zaga,df_72h_qrf,obs)


#72h_zaga_5

obs = as.data.frame(final2018_72h$R4final1)
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` <= 5)] <- 0
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` > 5)] <- 1

df_72h_zaga = cbind.data.frame(CDF_72h_zaga$V9,obs$`final2018_72h$R4final1`)

colnames(df_72h_zaga) = c("pred","obs")

#A<- verify(df_72h_zaga$obs, df_72h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 5mm threshold (72h forecast)")

brier_72h_zaga = brier(df_72h_zaga$obs,df_72h_zaga$pred, bins=F)
brier.ss_72h_zaga_5 = brier_72h_zaga$ss

#72h_qrf_5
obs = as.data.frame(final2018_72h$R4final1)
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` <= 5)] <- 0
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` > 5)] <- 1

df_72h_qrf = cbind.data.frame(CDF_72h_qrf$V9,obs$`final2018_72h$R4final1`)

colnames(df_72h_qrf) = c("pred","obs")
#A<- verify(df_72h_qrf$obs, df_72h_qrf$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of qrf with 5mm threshold (72h forecast)")


brier_72h_qrf = brier(df_72h_qrf$obs,df_72h_qrf$pred, bins=F)
brier.ss_72h_qrf_5 = brier_72h_qrf$ss
rm(brier_72h_zaga,brier_72h_qrf,df_72h_zaga,df_72h_qrf,obs)

#72h_zaga_10

obs = as.data.frame(final2018_72h$R4final1)
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` <= 10)] <- 0
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` >10)] <- 1

df_72h_zaga = cbind.data.frame(CDF_72h_zaga$V11,obs$`final2018_72h$R4final1`)

colnames(df_72h_zaga) = c("pred","obs")


#A<- verify(df_72h_zaga$obs, df_72h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 10mm threshold (72h forecast)")


brier_72h_zaga = brier(df_72h_zaga$obs,df_72h_zaga$pred, bins=F)
brier.ss_72h_zaga_10 = brier_72h_zaga$ss

#72h_qrf_10
obs = as.data.frame(final2018_72h$R4final1)
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` <=10)] <- 0
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` > 10)] <- 1

df_72h_qrf = cbind.data.frame(CDF_72h_qrf$V11,obs$`final2018_72h$R4final1`)

colnames(df_72h_qrf) = c("pred","obs")

#A<- verify(df_72h_qrf$obs, df_72h_qrf$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of qrf with 10mm threshold (72h forecast)")


brier_72h_qrf = brier(df_72h_qrf$obs,df_72h_qrf$pred, bins=F)
brier.ss_72h_qrf_10 = brier_72h_qrf$ss
rm(brier_72h_zaga,brier_72h_qrf,df_72h_zaga,df_72h_qrf,obs)



#72h_zaga_15

obs = as.data.frame(final2018_72h$R4final1)
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` <= 15)] <- 0
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` >15)] <- 1

df_72h_zaga = cbind.data.frame(CDF_72h_zaga$V13,obs$`final2018_72h$R4final1`)

colnames(df_72h_zaga) = c("pred","obs")


#A<- verify(df_72h_zaga$obs, df_72h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 15mm threshold (72h forecast)")

brier_72h_zaga = brier(df_72h_zaga$obs,df_72h_zaga$pred, bins=F)
brier.ss_72h_zaga_15 = brier_72h_zaga$ss

#72h_qrf_15
obs = as.data.frame(final2018_72h$R4final1)
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` <=15)] <- 0
obs$`final2018_72h$R4final1`[which(obs$`final2018_72h$R4final1` > 15)] <- 1

df_72h_qrf = cbind.data.frame(CDF_72h_qrf$V13,obs$`final2018_72h$R4final1`)

colnames(df_72h_qrf) = c("pred","obs")

#A<- verify(df_72h_qrf$obs, df_72h_qrf$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of qrf with 15mm threshold (72h forecast)")

brier_72h_qrf = brier(df_72h_qrf$obs,df_72h_qrf$pred, bins=F)
brier.ss_72h_qrf_15 = brier_72h_qrf$ss
rm(brier_72h_zaga,brier_72h_qrf,df_72h_zaga,df_72h_qrf,obs)

rm(final_72h,final2018_72h,final2020_72h)




#96h_zaga_0.05
CDF_96h_zaga = 1- CDF_96h_zaga
final_96h = read.table("E:/output/validating_model/96hours/96hour_forecast_alldata.txt", sep = "\t", header = T)

final_96h = na.omit(final_96h)
final_96h$R4final1 = as.numeric(final_96h$R4final1)
final2018_96h = final_96h[1:164126,]
final2020_96h = setdiff(final_96h,final2018_96h)

obs = as.data.frame(final2018_96h$R4final1)
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` <= 0.05)] <- 0
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` > 0.05)] <- 1

df_96h_zaga = cbind.data.frame(CDF_96h_zaga$V1,obs$`final2018_96h$R4final1`)

colnames(df_96h_zaga) = c("pred","obs")


#A<- verify(df_96h_zaga$obs, df_96h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 0.05mm threshold (96h forecast)")

brier_96h_zaga = brier(df_96h_zaga$obs,df_96h_zaga$pred, bins=F)
brier.ss_96h_zaga_0.05 = brier_96h_zaga$ss

#96h_qrf_0.05
CDF_96h_qrf = 1- CDF_96h_qrf
obs = as.data.frame(final2018_96h$R4final1)
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` <= 0.05)] <- 0
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` > 0.05)] <- 1

df_96h_qrf = cbind.data.frame(CDF_96h_qrf$V1,obs$`final2018_96h$R4final1`)

colnames(df_96h_qrf) = c("pred","obs")

#A<- verify(df_96h_qrf$obs, df_96h_qrf$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of qrf with 0.05mm threshold (96h forecast)")

brier_96h_qrf = brier(df_96h_qrf$obs,df_96h_qrf$pred, bins=F)
brier.ss_96h_qrf_0.05 = brier_96h_qrf$ss
rm(brier_96h_zaga,brier_96h_qrf,df_96h_zaga,df_96h_qrf,obs)

#96h_zaga_0.5

obs = as.data.frame(final2018_96h$R4final1)
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` <= 0.5)] <- 0
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` > 0.5)] <- 1

df_96h_zaga = cbind.data.frame(CDF_96h_zaga$V4,obs$`final2018_96h$R4final1`)

colnames(df_96h_zaga) = c("pred","obs")

brier_96h_zaga = brier(df_96h_zaga$obs,df_96h_zaga$pred, bins=F)
brier.ss_96h_zaga_0.5 = brier_96h_zaga$ss

#96h_qrf_0.5
obs = as.data.frame(final2018_96h$R4final1)
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` <= 0.5)] <- 0
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` > 0.5)] <- 1

df_96h_qrf = cbind.data.frame(CDF_96h_qrf$V4,obs$`final2018_96h$R4final1`)

colnames(df_96h_qrf) = c("pred","obs")

brier_96h_qrf = brier(df_96h_qrf$obs,df_96h_qrf$pred, bins=F)
brier.ss_96h_qrf_0.5 = brier_96h_qrf$ss
rm(brier_96h_zaga,brier_96h_qrf,df_96h_zaga,df_96h_qrf,obs)

#96h_zaga_1

obs = as.data.frame(final2018_96h$R4final1)
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` <= 1)] <- 0
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` > 1)] <- 1

df_96h_zaga = cbind.data.frame(CDF_96h_zaga$V6,obs$`final2018_96h$R4final1`)

colnames(df_96h_zaga) = c("pred","obs")

brier_96h_zaga = brier(df_96h_zaga$obs,df_96h_zaga$pred, bins=F)
brier.ss_96h_zaga_1 = brier_96h_zaga$ss

#96h_qrf_1
obs = as.data.frame(final2018_96h$R4final1)
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` <= 1)] <- 0
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` > 1)] <- 1

df_96h_qrf = cbind.data.frame(CDF_96h_qrf$V6,obs$`final2018_96h$R4final1`)

colnames(df_96h_qrf) = c("pred","obs")

brier_96h_qrf = brier(df_96h_qrf$obs,df_96h_qrf$pred, bins=F)
brier.ss_96h_qrf_1 = brier_96h_qrf$ss
rm(brier_96h_zaga,brier_96h_qrf,df_96h_zaga,df_96h_qrf,obs)


#96h_zaga_5

obs = as.data.frame(final2018_96h$R4final1)
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` <= 5)] <- 0
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` > 5)] <- 1

df_96h_zaga = cbind.data.frame(CDF_96h_zaga$V9,obs$`final2018_96h$R4final1`)

colnames(df_96h_zaga) = c("pred","obs")

#A<- verify(df_96h_zaga$obs, df_96h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 5mm threshold (96h forecast)")

brier_96h_zaga = brier(df_96h_zaga$obs,df_96h_zaga$pred, bins=F)
brier.ss_96h_zaga_5 = brier_96h_zaga$ss

#96h_qrf_5
obs = as.data.frame(final2018_96h$R4final1)
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` <= 5)] <- 0
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` > 5)] <- 1

df_96h_qrf = cbind.data.frame(CDF_96h_qrf$V9,obs$`final2018_96h$R4final1`)

colnames(df_96h_qrf) = c("pred","obs")
#A<- verify(df_96h_qrf$obs, df_96h_qrf$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of qrf with 5mm threshold (96h forecast)")


brier_96h_qrf = brier(df_96h_qrf$obs,df_96h_qrf$pred, bins=F)
brier.ss_96h_qrf_5 = brier_96h_qrf$ss
rm(brier_96h_zaga,brier_96h_qrf,df_96h_zaga,df_96h_qrf,obs)

#96h_zaga_10

obs = as.data.frame(final2018_96h$R4final1)
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` <= 10)] <- 0
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` >10)] <- 1

df_96h_zaga = cbind.data.frame(CDF_96h_zaga$V11,obs$`final2018_96h$R4final1`)

colnames(df_96h_zaga) = c("pred","obs")


#A<- verify(df_96h_zaga$obs, df_96h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 10mm threshold (96h forecast)")


brier_96h_zaga = brier(df_96h_zaga$obs,df_96h_zaga$pred, bins=F)
brier.ss_96h_zaga_10 = brier_96h_zaga$ss

#96h_qrf_10
obs = as.data.frame(final2018_96h$R4final1)
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` <=10)] <- 0
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` > 10)] <- 1

df_96h_qrf = cbind.data.frame(CDF_96h_qrf$V11,obs$`final2018_96h$R4final1`)

colnames(df_96h_qrf) = c("pred","obs")

#A<- verify(df_96h_qrf$obs, df_96h_qrf$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of qrf with 10mm threshold (96h forecast)")


brier_96h_qrf = brier(df_96h_qrf$obs,df_96h_qrf$pred, bins=F)
brier.ss_96h_qrf_10 = brier_96h_qrf$ss
rm(brier_96h_zaga,brier_96h_qrf,df_96h_zaga,df_96h_qrf,obs)



#96h_zaga_15

obs = as.data.frame(final2018_96h$R4final1)
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` <= 15)] <- 0
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` >15)] <- 1

df_96h_zaga = cbind.data.frame(CDF_96h_zaga$V13,obs$`final2018_96h$R4final1`)

colnames(df_96h_zaga) = c("pred","obs")


#A<- verify(df_96h_zaga$obs, df_96h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 15mm threshold (96h forecast)")

brier_96h_zaga = brier(df_96h_zaga$obs,df_96h_zaga$pred, bins=F)
brier.ss_96h_zaga_15 = brier_96h_zaga$ss

#96h_qrf_15
obs = as.data.frame(final2018_96h$R4final1)
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` <=15)] <- 0
obs$`final2018_96h$R4final1`[which(obs$`final2018_96h$R4final1` > 15)] <- 1

df_96h_qrf = cbind.data.frame(CDF_96h_qrf$V13,obs$`final2018_96h$R4final1`)

colnames(df_96h_qrf) = c("pred","obs")

#A<- verify(df_96h_qrf$obs, df_96h_qrf$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of qrf with 15mm threshold (96h forecast)")

brier_96h_qrf = brier(df_96h_qrf$obs,df_96h_qrf$pred, bins=F)
brier.ss_96h_qrf_15 = brier_96h_qrf$ss
rm(brier_96h_zaga,brier_96h_qrf,df_96h_zaga,df_96h_qrf,obs)

rm(final_96h,final2018_96h,final2020_96h)





#120h_zaga_0.05
CDF_120h_zaga = 1- CDF_120h_zaga
final_120h = read.table("E:/output/validating_model/120hours/120hour_forecast_alldata.txt", sep = "\t", header = T)

final_120h = na.omit(final_120h)
final_120h$R4final1 = as.numeric(final_120h$R4final1)
final2018_120h = final_120h[1:164126,]
final2020_120h = setdiff(final_120h,final2018_120h)

obs = as.data.frame(final2018_120h$R4final1)
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` <= 0.05)] <- 0
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` > 0.05)] <- 1

df_120h_zaga = cbind.data.frame(CDF_120h_zaga$V1,obs$`final2018_120h$R4final1`)

colnames(df_120h_zaga) = c("pred","obs")


#A<- verify(df_120h_zaga$obs, df_120h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 0.05mm threshold (120h forecast)")

brier_120h_zaga = brier(df_120h_zaga$obs,df_120h_zaga$pred, bins=F)
brier.ss_120h_zaga_0.05 = brier_120h_zaga$ss

#120h_qrf_0.05
CDF_120h_qrf = 1- CDF_120h_qrf
obs = as.data.frame(final2018_120h$R4final1)
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` <= 0.05)] <- 0
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` > 0.05)] <- 1

df_120h_qrf = cbind.data.frame(CDF_120h_qrf$V1,obs$`final2018_120h$R4final1`)

colnames(df_120h_qrf) = c("pred","obs")

#A<- verify(df_120h_qrf$obs, df_120h_qrf$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of qrf with 0.05mm threshold (120h forecast)")

brier_120h_qrf = brier(df_120h_qrf$obs,df_120h_qrf$pred, bins=F)
brier.ss_120h_qrf_0.05 = brier_120h_qrf$ss
rm(brier_120h_zaga,brier_120h_qrf,df_120h_zaga,df_120h_qrf,obs)

#120h_zaga_0.5

obs = as.data.frame(final2018_120h$R4final1)
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` <= 0.5)] <- 0
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` > 0.5)] <- 1

df_120h_zaga = cbind.data.frame(CDF_120h_zaga$V4,obs$`final2018_120h$R4final1`)

colnames(df_120h_zaga) = c("pred","obs")

brier_120h_zaga = brier(df_120h_zaga$obs,df_120h_zaga$pred, bins=F)
brier.ss_120h_zaga_0.5 = brier_120h_zaga$ss

#120h_qrf_0.5
obs = as.data.frame(final2018_120h$R4final1)
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` <= 0.5)] <- 0
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` > 0.5)] <- 1

df_120h_qrf = cbind.data.frame(CDF_120h_qrf$V4,obs$`final2018_120h$R4final1`)

colnames(df_120h_qrf) = c("pred","obs")

brier_120h_qrf = brier(df_120h_qrf$obs,df_120h_qrf$pred, bins=F)
brier.ss_120h_qrf_0.5 = brier_120h_qrf$ss
rm(brier_120h_zaga,brier_120h_qrf,df_120h_zaga,df_120h_qrf,obs)

#120h_zaga_1

obs = as.data.frame(final2018_120h$R4final1)
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` <= 1)] <- 0
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` > 1)] <- 1

df_120h_zaga = cbind.data.frame(CDF_120h_zaga$V6,obs$`final2018_120h$R4final1`)

colnames(df_120h_zaga) = c("pred","obs")

brier_120h_zaga = brier(df_120h_zaga$obs,df_120h_zaga$pred, bins=F)
brier.ss_120h_zaga_1 = brier_120h_zaga$ss

#120h_qrf_1
obs = as.data.frame(final2018_120h$R4final1)
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` <= 1)] <- 0
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` > 1)] <- 1

df_120h_qrf = cbind.data.frame(CDF_120h_qrf$V6,obs$`final2018_120h$R4final1`)

colnames(df_120h_qrf) = c("pred","obs")

brier_120h_qrf = brier(df_120h_qrf$obs,df_120h_qrf$pred, bins=F)
brier.ss_120h_qrf_1 = brier_120h_qrf$ss
rm(brier_120h_zaga,brier_120h_qrf,df_120h_zaga,df_120h_qrf,obs)


#120h_zaga_5

obs = as.data.frame(final2018_120h$R4final1)
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` <= 5)] <- 0
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` > 5)] <- 1

df_120h_zaga = cbind.data.frame(CDF_120h_zaga$V9,obs$`final2018_120h$R4final1`)

colnames(df_120h_zaga) = c("pred","obs")

#A<- verify(df_120h_zaga$obs, df_120h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 5mm threshold (120h forecast)")

brier_120h_zaga = brier(df_120h_zaga$obs,df_120h_zaga$pred, bins=F)
brier.ss_120h_zaga_5 = brier_120h_zaga$ss

#120h_qrf_5
obs = as.data.frame(final2018_120h$R4final1)
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` <= 5)] <- 0
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` > 5)] <- 1

df_120h_qrf = cbind.data.frame(CDF_120h_qrf$V9,obs$`final2018_120h$R4final1`)

colnames(df_120h_qrf) = c("pred","obs")
#A<- verify(df_120h_qrf$obs, df_120h_qrf$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of qrf with 5mm threshold (120h forecast)")


brier_120h_qrf = brier(df_120h_qrf$obs,df_120h_qrf$pred, bins=F)
brier.ss_120h_qrf_5 = brier_120h_qrf$ss
rm(brier_120h_zaga,brier_120h_qrf,df_120h_zaga,df_120h_qrf,obs)

#120h_zaga_10

obs = as.data.frame(final2018_120h$R4final1)
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` <= 10)] <- 0
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` >10)] <- 1

df_120h_zaga = cbind.data.frame(CDF_120h_zaga$V11,obs$`final2018_120h$R4final1`)

colnames(df_120h_zaga) = c("pred","obs")


#A<- verify(df_120h_zaga$obs, df_120h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 10mm threshold (120h forecast)")


brier_120h_zaga = brier(df_120h_zaga$obs,df_120h_zaga$pred, bins=F)
brier.ss_120h_zaga_10 = brier_120h_zaga$ss

#120h_qrf_10
obs = as.data.frame(final2018_120h$R4final1)
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` <=10)] <- 0
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` > 10)] <- 1

df_120h_qrf = cbind.data.frame(CDF_120h_qrf$V11,obs$`final2018_120h$R4final1`)

colnames(df_120h_qrf) = c("pred","obs")

#A<- verify(df_120h_qrf$obs, df_120h_qrf$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of qrf with 10mm threshold (120h forecast)")


brier_120h_qrf = brier(df_120h_qrf$obs,df_120h_qrf$pred, bins=F)
brier.ss_120h_qrf_10 = brier_120h_qrf$ss
rm(brier_120h_zaga,brier_120h_qrf,df_120h_zaga,df_120h_qrf,obs)



#120h_zaga_15

obs = as.data.frame(final2018_120h$R4final1)
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` <= 15)] <- 0
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` >15)] <- 1

df_120h_zaga = cbind.data.frame(CDF_120h_zaga$V13,obs$`final2018_120h$R4final1`)

colnames(df_120h_zaga) = c("pred","obs")


#A<- verify(df_120h_zaga$obs, df_120h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 15mm threshold (120h forecast)")

brier_120h_zaga = brier(df_120h_zaga$obs,df_120h_zaga$pred, bins=F)
brier.ss_120h_zaga_15 = brier_120h_zaga$ss

#120h_qrf_15
obs = as.data.frame(final2018_120h$R4final1)
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` <=15)] <- 0
obs$`final2018_120h$R4final1`[which(obs$`final2018_120h$R4final1` > 15)] <- 1

df_120h_qrf = cbind.data.frame(CDF_120h_qrf$V13,obs$`final2018_120h$R4final1`)

colnames(df_120h_qrf) = c("pred","obs")

#A<- verify(df_120h_qrf$obs, df_120h_qrf$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of qrf with 15mm threshold (120h forecast)")

brier_120h_qrf = brier(df_120h_qrf$obs,df_120h_qrf$pred, bins=F)
brier.ss_120h_qrf_15 = brier_120h_qrf$ss
rm(brier_120h_zaga,brier_120h_qrf,df_120h_zaga,df_120h_qrf,obs)

rm(final_120h,final2018_120h,final2020_120h)






#144h_zaga_0.05
CDF_144h_zaga = 1- CDF_144h_zaga
final_144h = read.table("E:/output/validating_model/144hours/144hour_forecast_alldata.txt", sep = "\t", header = T)

final_144h = na.omit(final_144h)
final_144h$R4final1 = as.numeric(final_144h$R4final1)
final2018_144h = final_144h[1:164126,]
final2020_144h = setdiff(final_144h,final2018_144h)

obs = as.data.frame(final2018_144h$R4final1)
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` <= 0.05)] <- 0
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` > 0.05)] <- 1

df_144h_zaga = cbind.data.frame(CDF_144h_zaga$V1,obs$`final2018_144h$R4final1`)

colnames(df_144h_zaga) = c("pred","obs")


A<- verify(df_144h_zaga$obs, df_144h_zaga$pred, frcst.type = "prob", obs.type = "binary")

reliability.plot(A, titl = "Reliability diagram of zaga with 0.05mm threshold (144h forecast)")

brier_144h_zaga = brier(df_144h_zaga$obs,df_144h_zaga$pred, bins=F)
brier.ss_144h_zaga_0.05 = brier_144h_zaga$ss

#144h_qrf_0.05
CDF_144h_qrf = 1- CDF_144h_qrf
obs = as.data.frame(final2018_144h$R4final1)
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` <= 0.05)] <- 0
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` > 0.05)] <- 1

df_144h_qrf = cbind.data.frame(CDF_144h_qrf$V1,obs$`final2018_144h$R4final1`)

colnames(df_144h_qrf) = c("pred","obs")

A<- verify(df_144h_qrf$obs, df_144h_qrf$pred, frcst.type = "prob", obs.type = "binary")

reliability.plot(A, titl = "Reliability diagram of qrf with 0.05mm threshold (144h forecast)")

brier_144h_qrf = brier(df_144h_qrf$obs,df_144h_qrf$pred, bins=F)
brier.ss_144h_qrf_0.05 = brier_144h_qrf$ss
rm(brier_144h_zaga,brier_144h_qrf,df_144h_zaga,df_144h_qrf,obs)

#144h_zaga_0.5

obs = as.data.frame(final2018_144h$R4final1)
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` <= 0.5)] <- 0
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` > 0.5)] <- 1

df_144h_zaga = cbind.data.frame(CDF_144h_zaga$V4,obs$`final2018_144h$R4final1`)

colnames(df_144h_zaga) = c("pred","obs")

A<- verify(df_144h_zaga$obs, df_144h_zaga$pred, frcst.type = "prob", obs.type = "binary")

reliability.plot(A, titl = "Reliability diagram of zaga with 1mm threshold (144h forecast)")


brier_144h_zaga = brier(df_144h_zaga$obs,df_144h_zaga$pred, bins=F)
brier.ss_144h_zaga_0.5 = brier_144h_zaga$ss

#144h_qrf_0.5
obs = as.data.frame(final2018_144h$R4final1)
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` <= 0.5)] <- 0
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` > 0.5)] <- 1

df_144h_qrf = cbind.data.frame(CDF_144h_qrf$V4,obs$`final2018_144h$R4final1`)

colnames(df_144h_qrf) = c("pred","obs")

A<- verify(df_144h_qrf$obs, df_144h_qrf$pred, frcst.type = "prob", obs.type = "binary")

reliability.plot(A, titl = "Reliability diagram of qrf with 1mm threshold (144h forecast)")


brier_144h_qrf = brier(df_144h_qrf$obs,df_144h_qrf$pred, bins=F)
brier.ss_144h_qrf_0.5 = brier_144h_qrf$ss
rm(brier_144h_zaga,brier_144h_qrf,df_144h_zaga,df_144h_qrf,obs)

#144h_zaga_1

obs = as.data.frame(final2018_144h$R4final1)
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` <= 1)] <- 0
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` > 1)] <- 1

df_144h_zaga = cbind.data.frame(CDF_144h_zaga$V6,obs$`final2018_144h$R4final1`)

colnames(df_144h_zaga) = c("pred","obs")


A<- verify(df_144h_zaga$obs, df_144h_zaga$pred, frcst.type = "prob", obs.type = "binary")

reliability.plot(A, titl = "Reliability diagram of zaga with 1mm threshold (144h forecast)")

brier_144h_zaga = brier(df_144h_zaga$obs,df_144h_zaga$pred, bins=F)
brier.ss_144h_zaga_1 = brier_144h_zaga$ss

#144h_qrf_1
obs = as.data.frame(final2018_144h$R4final1)
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` <= 1)] <- 0
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` > 1)] <- 1

df_144h_qrf = cbind.data.frame(CDF_144h_qrf$V6,obs$`final2018_144h$R4final1`)

colnames(df_144h_qrf) = c("pred","obs")

A<- verify(df_144h_qrf$obs, df_144h_qrf$pred, frcst.type = "prob", obs.type = "binary")

reliability.plot(A, titl = "Reliability diagram of qrf with 1mm threshold (144h forecast)")


brier_144h_qrf = brier(df_144h_qrf$obs,df_144h_qrf$pred, bins=F)
brier.ss_144h_qrf_1 = brier_144h_qrf$ss
rm(brier_144h_zaga,brier_144h_qrf,df_144h_zaga,df_144h_qrf,obs)


#144h_zaga_5

obs = as.data.frame(final2018_144h$R4final1)
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` <= 5)] <- 0
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` > 5)] <- 1

df_144h_zaga = cbind.data.frame(CDF_144h_zaga$V9,obs$`final2018_144h$R4final1`)

colnames(df_144h_zaga) = c("pred","obs")

#A<- verify(df_144h_zaga$obs, df_144h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 5mm threshold (144h forecast)")

brier_144h_zaga = brier(df_144h_zaga$obs,df_144h_zaga$pred, bins=F)
brier.ss_144h_zaga_5 = brier_144h_zaga$ss

#144h_qrf_5
obs = as.data.frame(final2018_144h$R4final1)
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` <= 5)] <- 0
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` > 5)] <- 1

df_144h_qrf = cbind.data.frame(CDF_144h_qrf$V9,obs$`final2018_144h$R4final1`)

colnames(df_144h_qrf) = c("pred","obs")
#A<- verify(df_144h_qrf$obs, df_144h_qrf$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of qrf with 5mm threshold (144h forecast)")


brier_144h_qrf = brier(df_144h_qrf$obs,df_144h_qrf$pred, bins=F)
brier.ss_144h_qrf_5 = brier_144h_qrf$ss
rm(brier_144h_zaga,brier_144h_qrf,df_144h_zaga,df_144h_qrf,obs)

#144h_zaga_10

obs = as.data.frame(final2018_144h$R4final1)
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` <= 10)] <- 0
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` >10)] <- 1

df_144h_zaga = cbind.data.frame(CDF_144h_zaga$V11,obs$`final2018_144h$R4final1`)

colnames(df_144h_zaga) = c("pred","obs")


#A<- verify(df_144h_zaga$obs, df_144h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 10mm threshold (144h forecast)")


brier_144h_zaga = brier(df_144h_zaga$obs,df_144h_zaga$pred, bins=F)
brier.ss_144h_zaga_10 = brier_144h_zaga$ss

#144h_qrf_10
obs = as.data.frame(final2018_144h$R4final1)
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` <=10)] <- 0
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` > 10)] <- 1

df_144h_qrf = cbind.data.frame(CDF_144h_qrf$V11,obs$`final2018_144h$R4final1`)

colnames(df_144h_qrf) = c("pred","obs")

#A<- verify(df_144h_qrf$obs, df_144h_qrf$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of qrf with 10mm threshold (144h forecast)")


brier_144h_qrf = brier(df_144h_qrf$obs,df_144h_qrf$pred, bins=F)
brier.ss_144h_qrf_10 = brier_144h_qrf$ss
rm(brier_144h_zaga,brier_144h_qrf,df_144h_zaga,df_144h_qrf,obs)



#144h_zaga_15

obs = as.data.frame(final2018_144h$R4final1)
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` <= 15)] <- 0
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` >15)] <- 1

df_144h_zaga = cbind.data.frame(CDF_144h_zaga$V13,obs$`final2018_144h$R4final1`)

colnames(df_144h_zaga) = c("pred","obs")


#A<- verify(df_144h_zaga$obs, df_144h_zaga$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of zaga with 15mm threshold (144h forecast)")

brier_144h_zaga = brier(df_144h_zaga$obs,df_144h_zaga$pred, bins=F)
brier.ss_144h_zaga_15 = brier_144h_zaga$ss

#144h_qrf_15
obs = as.data.frame(final2018_144h$R4final1)
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` <=15)] <- 0
obs$`final2018_144h$R4final1`[which(obs$`final2018_144h$R4final1` > 15)] <- 1

df_144h_qrf = cbind.data.frame(CDF_144h_qrf$V13,obs$`final2018_144h$R4final1`)

colnames(df_144h_qrf) = c("pred","obs")

#A<- verify(df_144h_qrf$obs, df_144h_qrf$pred, frcst.type = "prob", obs.type = "binary")

#reliability.plot(A, titl = "Reliability diagram of qrf with 15mm threshold (144h forecast)")

brier_144h_qrf = brier(df_144h_qrf$obs,df_144h_qrf$pred, bins=F)
brier.ss_144h_qrf_15 = brier_144h_qrf$ss
rm(brier_144h_zaga,brier_144h_qrf,df_144h_zaga,df_144h_qrf,obs)

rm(final_144h,final2018_144h,final2020_144h)

#df for 0.05 threshold
x <- rep(c(24,48,72,96,120,144), each = 2)
y <- rep(c('zaga','qrf'),times = 6)
value1 <- c(brier.ss_24h_zaga_0.05,brier.ss_24h_qrf_0.05,brier.ss_48h_zaga_0.05,brier.ss_48h_qrf_0.05,
            brier.ss_72h_zaga_0.05,brier.ss_72h_qrf_0.05,brier.ss_96h_zaga_0.05,brier.ss_96h_qrf_0.05,
            brier.ss_120h_zaga_0.05,brier.ss_120h_qrf_0.05,brier.ss_144h_zaga_0.05,brier.ss_144h_qrf_0.05)



#df for 0.5 threshold

value2 <- c(brier.ss_24h_zaga_0.5,brier.ss_24h_qrf_0.5,brier.ss_48h_zaga_0.5,brier.ss_48h_qrf_0.5,
            brier.ss_72h_zaga_0.5,brier.ss_72h_qrf_0.5,brier.ss_96h_zaga_0.5,brier.ss_96h_qrf_0.5,
            brier.ss_120h_zaga_0.5,brier.ss_120h_qrf_0.5,brier.ss_144h_zaga_0.5,brier.ss_144h_qrf_0.5)


#df for 1 threshold
value3 <- c(brier.ss_24h_zaga_1,brier.ss_24h_qrf_1,brier.ss_48h_zaga_1,brier.ss_48h_qrf_1,
            brier.ss_72h_zaga_1,brier.ss_72h_qrf_1,brier.ss_96h_zaga_1,brier.ss_96h_qrf_1,
            brier.ss_120h_zaga_1,brier.ss_120h_qrf_1,brier.ss_144h_zaga_1,brier.ss_144h_qrf_1)


#df for 5 threshold
value4 <- c(brier.ss_24h_zaga_5,brier.ss_24h_qrf_5,brier.ss_48h_zaga_5,brier.ss_48h_qrf_5,
            brier.ss_72h_zaga_5,brier.ss_72h_qrf_5,brier.ss_96h_zaga_5,brier.ss_96h_qrf_5,
            brier.ss_120h_zaga_5,brier.ss_120h_qrf_5,brier.ss_144h_zaga_5,brier.ss_144h_qrf_5)


#df for 10 threshold
value5 <- c(brier.ss_24h_zaga_10,brier.ss_24h_qrf_10,brier.ss_48h_zaga_10,brier.ss_48h_qrf_10,
            brier.ss_72h_zaga_10,brier.ss_72h_qrf_10,brier.ss_96h_zaga_10,brier.ss_96h_qrf_10,
            brier.ss_120h_zaga_10,brier.ss_120h_qrf_10,brier.ss_144h_zaga_10,brier.ss_144h_qrf_10)



value6 <- c(brier.ss_24h_zaga_15,brier.ss_24h_qrf_15,brier.ss_48h_zaga_15,brier.ss_48h_qrf_15,
            brier.ss_72h_zaga_15,brier.ss_72h_qrf_15,brier.ss_96h_zaga_15,brier.ss_96h_qrf_15,
            brier.ss_120h_zaga_15,brier.ss_120h_qrf_15,brier.ss_144h_zaga_15,brier.ss_144h_qrf_15)

df <- data.frame(x = x, y=y, value1=value1,value2=value2,value3=value3,value4=value4,
                 value5=value5,value6=value6)



p1<-ggplot(data = df, mapping = aes(x = factor(x), y = value1, fill = y)) + geom_bar(stat = 'identity', position = 'dodge')+xlab("lead time (h)")+ylab("brier skill score")+
  ggtitle("threshold: 0.05mm/6h") +
  scale_fill_discrete(name = "Method")+
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))

p2<-ggplot(data = df, mapping = aes(x = factor(x), y = value2, fill = y)) + geom_bar(stat = 'identity', position = 'dodge')+xlab("lead time (h)")+ylab("brier skill score")+
  ggtitle("threshold: 0.5mm/6h") +
  scale_fill_discrete(name = "Method")+
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))


p3<-ggplot(data = df, mapping = aes(x = factor(x), y = value3, fill = y)) + geom_bar(stat = 'identity', position = 'dodge')+xlab("lead time (h)")+ylab("brier skill score")+
  ggtitle("threshold: 1mm/6h") +  
  scale_fill_discrete(name = "Method")+
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))

p4<-ggplot(data = df, mapping = aes(x = factor(x), y = value4, fill = y)) + geom_bar(stat = 'identity', position = 'dodge')+xlab("lead time (h)")+ylab("brier skill score")+
  ggtitle("threshold: 5mm/6h") +
  scale_fill_discrete(name = "Method")+
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))

p5<-ggplot(data = df, mapping = aes(x = factor(x), y = value5, fill = y)) + geom_bar(stat = 'identity', position = 'dodge')+xlab("lead time (h)")+ylab("brier skill score")+
  ggtitle("threshold: 10mm/6h") +
  scale_fill_discrete(name = "Method")+
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))

p6<-ggplot(data = df, mapping = aes(x = factor(x), y = value6, fill = y)) + geom_bar(stat = 'identity', position = 'dodge')+xlab("lead time (h)")+ylab("brier skill score")+
  ggtitle("threshold: 15mm/6h") +
  scale_fill_discrete(name = "Method")+
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))

ggarrange(p1,p2,p3,p4,p5,p6,ncol = 2,nrow=3,common.legend = T,legend = "right")

