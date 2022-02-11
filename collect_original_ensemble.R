######Collect original ensemble data#####################################################################
library(e1071)
Sys.time()
setwd('E:/ECMWFdata/ECME_2019')
first_category_name = list.files("Jan")
dir = paste("./Jan/",first_category_name,sep="")
n = length(dir) 
n_sub<-rep(0,n)
n_sub<-as.data.frame(n_sub)
n_sub<-t(n_sub)
head(n_sub)
ensemble_24hours = as.data.frame(matrix(NA, ncol = 51))


for(i in 1:n){         #对于每个一级目录(文件夹)
  b=list.files(dir[i]) #b是列出每个一级目录(文件夹)中每个xlsx文件的名称
  n_sub[i]=length(b) 

for(j in 1:n_sub[i]){     #对于每个一级目录(文件夹)下的每个xlsx文件
  data<-read.delim(file=paste(dir[i],'/',b[j],sep=''),sep = "", header = F) #读取xlsx文件
  #RRcum (mm)
  data_RRcum = data[479:531,]
  data_RRcum[,48:83]=trunc(data_RRcum[,48:83]/1000000)
  data_RRcum[2,64]=data_RRcum[2,64]/10^234
  data_RRcum = as.data.frame(lapply(data_RRcum,as.numeric))
  data_RRcum[2:53,]=data_RRcum[2:53,]/10
  RRcum_6hstep = data_RRcum[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
  #RR
  RR_6hstep = RRcum_6hstep[,1:(ncol(RRcum_6hstep)-1)]
  RR_6hstep = cbind(0,RR_6hstep)
  RRcum_6hstep = as.data.frame(lapply(RRcum_6hstep,as.numeric))
  RR_6hstep = as.data.frame(lapply(RR_6hstep,as.numeric))
  RR_6hstep = RRcum_6hstep-RR_6hstep
  RR_ENS = RR_6hstep[3:53,]
  RR_ENS = as.data.frame(lapply(RR_ENS,as.numeric))
  RR = RR_ENS[,4]
  ensemble_24hours = rbind.data.frame(ensemble_24hours,RR)
  }    
}  
Sys.time()    
ensemble_24hours = ensemble_24hours[-1,]
ensemble_24hours_1 = rbind.data.frame(ensemble_24hours_1,ensemble_24hours)
ensemble_24hours_1 = read.table("E:/output/ensemble_24hours.txt", sep = "\t", header = T)
########Match the original ensemble data with the observation#################################
library(lubridate)
final2018 = final2018[,99:102]
 colnames(final2018) = c("final","dftime","grid","rID")
 ensemble_24hours_1$grid = as.numeric(substr(ensemble_24hours_1$second_category,12,15))
 ensemble_24hours_1$dftime = as.POSIXct("2000-1-1 01:00:00", tz = "UTC")
 year(ensemble_24hours_1$dftime) = as.numeric(substr(ensemble_24hours_1$first_category,10,13))
 month(ensemble_24hours_1$dftime) = as.numeric(substr(ensemble_24hours_1$first_category,14,15))
 mday(ensemble_24hours_1$dftime) = as.numeric(substr(ensemble_24hours_1$first_category,16,17))
 hour(ensemble_24hours_1$dftime) = as.numeric(substr(ensemble_24hours_1$first_category,18,19))
 minute(ensemble_24hours_1$dftime) = 00
 second(ensemble_24hours_1$dftime) = 00
 ensemble_24hours_1$dftime = as.character(ensemble_24hours_1$dftime)
check = left_join(final2018,ensemble_24hours_1,by = c("dftime","grid"))
check = check[,-57]
check = check[,-56]
######Calculate the CRPS###############################################
crps_1 = as.data.frame(matrix(NA,nrow = dim(final2018)[1],ncol = 1))
library(scoringRules)
for (i in 1:dim(final2018)[1]){
  sample = as.numeric(check[i,5:55])
  crps_1[i,]=crps_sample(y = check$final[i],dat = sample)
}