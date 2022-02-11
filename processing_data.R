

##########Preparation part########################################################################
new = read.table("E:/output/time_information.txt", sep = "\t", header = T)
#new$dates = with_tz(new$dates,tzone = "UTC")
#new$dates = new$dates +hours(1)
#Rdata = as.data.frame(t(Rdata))
Rdata = cbind.data.frame(new,t(Rdata))

library(lubridate)                             
upscal_match = read.table("E:/KNMI_DATA/upscal_match_9km.txt", sep = "", header = T)
ll = c("time","time")
mm = c("match","match")
upscal_match = rbind.data.frame(mm,ll,upscal_match)
colnames(Rdata) = c("Rtime",1:599)




RmatchID = as.data.frame(c(1:2424))
colnames(RmatchID) = c("RmatchID")
Rdata = cbind.data.frame(RmatchID,Rdata)
#colnames(Rdata) = c(1:dim(Rdata)[2])
#Match data temporally
#######Match EC-files and Radar files###############################################

library(tidyverse)
library(dplyr)



  PP = function(data){
    data$grid = as.numeric(substr(data$second_category,12,15))
    
    data$dftime = as.POSIXct("2000-1-1 01:00:00", tz = "UTC")
    year(data$dftime) = as.numeric(substr(data$first_category,10,13))
    month(data$dftime) = as.numeric(substr(data$first_category,14,15))
    mday(data$dftime) = as.numeric(substr(data$first_category,16,17))
    hour(data$dftime) = as.numeric(substr(data$first_category,18,19))
    minute(data$dftime) = 00
    second(data$dftime) = 00
    data$dftime = as.character(data$dftime)
    
    ECsub = data[,2:5]
    Rsub = Rdata[,1:2]
    ECmatch = left_join(ECsub,Rsub,by = c("dftime"="Rtime"))
    data = cbind.data.frame(data,ECmatch$RmatchID)
    data = data[which(data$`ECmatch$RmatchID` <665 | data$`ECmatch$RmatchID` >783 & data$`ECmatch$RmatchID`<1517
                      | data$`ECmatch$RmatchID`>1635 & data$`ECmatch$RmatchID`<2365 ),]   
    

  }
  
  result <- list(RRmean, RRsd, RRmedian, RRq10, RRq25, RRq75, RRq90, RRp0.05, RRp0.3, RRp10, RRp20, RRmin, RRmax, 
                 RRHRES, WSPDmean, WSPDsd, WSPDmedian,WSPDq10, WSPDq25, WSPDq75, WSPDq90, WSPDmin, WSPDmax, WSPDHRES, 
                 T2Mmean, T2Msd, T2Mmedian, T2Mq10, T2Mq25, T2Mq75, T2Mq90, T2Mmin, T2Mmax, T2MHRES, EVAmean, EVAsd, 
                 EVAmedian, EVAq10, EVAq25, EVAq75, EVAq90, EVAmin, EVAmax, EVAHRES, DPTmean, DPTmedian, DPTsd, DPTq10,
                 DPTq25, DPTq75, DPTq90, DPTmin, DPTmax, DPTHRES, CONRRmean, CONRRmedian, CONRRsd, CONRRq10, CONRRq25, 
                 CONRRq75, CONRRq90, CONRRp0.05, CONRRp0.3, CONRRp10, CONRRp20, CONRRmin, CONRRmax, CONRRHRES, CCmean, 
                 CCmedian, CCsd, CCq10, CCq25, CCq75, CCq90, CCmin, CCmax, CCHRES, CAPEmean, CAPEmedian, CAPEsd, 
                 CAPEq10, CAPEq25, CAPEq75, CAPEq90, CAPEmin, CAPEmax, CAPEHRES, CAPESmean, CAPESmedian, CAPESsd, 
                 CAPESq10, CAPESq25, CAPESq75, CAPESq90, CAPESmin, CAPESmax, CAPESHRES) %>%  lapply(PP)  
   
  
  RRmean = as.data.frame(result[1])
  RRsd =  as.data.frame(result[2])
  RRmedian =  as.data.frame(result[3]) 
  RRq10 = as.data.frame(result[4])
  RRq25 =  as.data.frame(result[5])
  RRq75 = as.data.frame(result[6])
  RRq90 = as.data.frame(result[7])
  RRp0.05 = as.data.frame(result[8])
  RRp0.3 = as.data.frame(result[9])
  RRp10 = as.data.frame(result[10])
  RRp20 = as.data.frame(result[11])
  RRmin = as.data.frame(result[12]) 
  RRmax = as.data.frame(result[13])
  RRHRES = as.data.frame(result[14])
  WSPDmean = as.data.frame(result[15])
  WSPDsd = as.data.frame(result[16])
  WSPDmedian = as.data.frame(result[17])
  WSPDq10 = as.data.frame(result[18]) 
  WSPDq25 = as.data.frame(result[19])
  WSPDq75 = as.data.frame(result[20])
  WSPDq90 = as.data.frame(result[21])
  WSPDmin = as.data.frame(result[22])
  WSPDmax = as.data.frame(result[23])
  WSPDHRES = as.data.frame(result[24]) 
  T2Mmean = as.data.frame(result[25])
  T2Msd = as.data.frame(result[26])
  T2Mmedian = as.data.frame(result[27])
  T2Mq10 = as.data.frame(result[28])
  T2Mq25 = as.data.frame(result[29])
  T2Mq75 = as.data.frame(result[30])
  T2Mq90 = as.data.frame(result[31])
  T2Mmin = as.data.frame(result[32])
  T2Mmax = as.data.frame(result[33])
  T2MHRES = as.data.frame(result[34])
  EVAmean = as.data.frame(result[35])
  EVAsd = as.data.frame(result[36])
  EVAmedian = as.data.frame(result[37])
  EVAq10 = as.data.frame(result[38])
  EVAq25 = as.data.frame(result[39])
  EVAq75 = as.data.frame(result[40])
  EVAq90 = as.data.frame(result[41])
  EVAmin = as.data.frame(result[42])
  EVAmax = as.data.frame(result[43])
  EVAHRES = as.data.frame(result[44])
  DPTmean = as.data.frame(result[45])
  DPTmedian = as.data.frame(result[46])
  DPTsd = as.data.frame(result[47])
  DPTq10 = as.data.frame(result[48])
  DPTq25 = as.data.frame(result[49])
  DPTq75 = as.data.frame(result[50])
  DPTq90 = as.data.frame(result[51])
  DPTmin = as.data.frame(result[52])
  DPTmax = as.data.frame(result[53])
  DPTHRES = as.data.frame(result[54])
  CONRRmean = as.data.frame(result[55])
  CONRRmedian = as.data.frame(result[56])
  CONRRsd = as.data.frame(result[57])
  CONRRq10 = as.data.frame(result[58])
  CONRRq25 = as.data.frame(result[59])
  CONRRq75 = as.data.frame(result[60])
  CONRRq90 = as.data.frame(result[61])
  CONRRp0.05 = as.data.frame(result[62])
  CONRRp0.3 = as.data.frame(result[63])
  CONRRp10 = as.data.frame(result[64])
  CONRRp20 = as.data.frame(result[65])
  CONRRmin = as.data.frame(result[66])
  CONRRmax = as.data.frame(result[67])
  CONRRHRES = as.data.frame(result[68])
  CCmean = as.data.frame(result[69])
  CCmedian = as.data.frame(result[70])
  CCsd = as.data.frame(result[71])
  CCq10 = as.data.frame(result[72])
  CCq25 = as.data.frame(result[73])
  CCq75 = as.data.frame(result[74])
  CCq90 = as.data.frame(result[75])
  CCmin = as.data.frame(result[76])
  CCmax = as.data.frame(result[77])
  CCHRES = as.data.frame(result[78])
  CAPEmean = as.data.frame(result[79])
  CAPEmedian = as.data.frame(result[80])
  CAPEsd = as.data.frame(result[81]) 
  CAPEq10 = as.data.frame(result[82])
  CAPEq25 = as.data.frame(result[83])
  CAPEq75 = as.data.frame(result[84])
  CAPEq90 = as.data.frame(result[85])
  CAPEmin = as.data.frame(result[86])
  CAPEmax = as.data.frame(result[87])
  CAPEHRES = as.data.frame(result[88])
  CAPESmean = as.data.frame(result[89])
  CAPESmedian = as.data.frame(result[90])
  CAPESsd = as.data.frame(result[91]) 
  CAPESq10 = as.data.frame(result[92])
  CAPESq25 = as.data.frame(result[93])
  CAPESq75 = as.data.frame(result[94])
  CAPESq90 = as.data.frame(result[95])
  CAPESmin = as.data.frame(result[96])
  CAPESmax = as.data.frame(result[97])
  CAPESHRES = as.data.frame(result[98])
  
 
  
   
  
  




RRmeanfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
RRmeanfinal <- RRmeanfinal[-(1:2),]
EC4time <- as.data.frame(matrix(upscal_match$rID,nrow=601,ncol=1))
EC4time <- EC4time[-(1:2),]
RRsdfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
RRsdfinal <- RRsdfinal[-(1:2),]
RRmedianfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
RRmedianfinal <- RRmedianfinal[-(1:2),]
RRq10final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
RRq10final <- RRq10final[-(1:2),]
RRq25final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
RRq25final <- RRq25final[-(1:2),]
RRq75final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
RRq75final <- RRq75final[-(1:2),]
RRq90final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
RRq90final <- RRq90final[-(1:2),]
RRP0.05final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
RRP0.05final <- RRP0.05final[-(1:2),]
RRP0.3final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
RRP0.3final <- RRP0.3final[-(1:2),]
RRP10final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
RRP10final <- RRP10final[-(1:2),]
RRP20final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
RRP20final <- RRP20final[-(1:2),]
RRminfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
RRminfinal <- RRminfinal[-(1:2),]
RRmaxfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
RRmaxfinal <- RRmaxfinal[-(1:2),]
RRHRESfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
RRHRESfinal <- RRHRESfinal[-(1:2),]
R4final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
R4final <- R4final[-(1:2),]
EC4ID <- as.data.frame(matrix(upscal_match$rID,nrow = 601,ncol=1))
EC4ID <- EC4ID[-(1:2),]
R4ID <- as.data.frame(matrix(upscal_match$rID,nrow = 601,ncol=1))
R4ID <- R4ID[-(1:2),]


CONRRmeanfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CONRRmeanfinal <- CONRRmeanfinal[-(1:2),]
CONRRsdfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CONRRsdfinal <- CONRRsdfinal[-(1:2),]
CONRRmedianfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CONRRmedianfinal <- CONRRmedianfinal[-(1:2),]
CONRRq10final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CONRRq10final <- CONRRq10final[-(1:2),]
CONRRq25final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CONRRq25final <- CONRRq25final[-(1:2),]
CONRRq75final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CONRRq75final <- CONRRq75final[-(1:2),]
CONRRq90final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CONRRq90final <- CONRRq90final[-(1:2),]
CONRRP0.05final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CONRRP0.05final <- CONRRP0.05final[-(1:2),]
CONRRP0.3final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CONRRP0.3final <- CONRRP0.3final[-(1:2),]
CONRRP10final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CONRRP10final <- CONRRP10final[-(1:2),]
CONRRP20final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CONRRP20final <- CONRRP20final[-(1:2),]
CONRRminfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CONRRminfinal <- CONRRminfinal[-(1:2),]
CONRRmaxfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CONRRmaxfinal <- CONRRmaxfinal[-(1:2),]
CONRRHRESfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CONRRHRESfinal <- CONRRHRESfinal[-(1:2),]

WSPDmeanfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
WSPDmeanfinal <- WSPDmeanfinal[-(1:2),]
WSPDsdfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
WSPDsdfinal <- WSPDsdfinal[-(1:2),]
WSPDmedianfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
WSPDmedianfinal <- WSPDmedianfinal[-(1:2),]
WSPDq10final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
WSPDq10final <- WSPDq10final[-(1:2),]
WSPDq25final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
WSPDq25final <- WSPDq25final[-(1:2),]
WSPDq75final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
WSPDq75final <- WSPDq75final[-(1:2),]
WSPDq90final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
WSPDq90final <- WSPDq90final[-(1:2),]
WSPDminfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
WSPDminfinal <- WSPDminfinal[-(1:2),]
WSPDmaxfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
WSPDmaxfinal <- WSPDmaxfinal[-(1:2),]
WSPDHRESfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
WSPDHRESfinal <- WSPDHRESfinal[-(1:2),]


DPTmeanfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
DPTmeanfinal <- DPTmeanfinal[-(1:2),]
DPTsdfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
DPTsdfinal <- DPTsdfinal[-(1:2),]
DPTmedianfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
DPTmedianfinal <- DPTmedianfinal[-(1:2),]
DPTq10final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
DPTq10final <- DPTq10final[-(1:2),]
DPTq25final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
DPTq25final <- DPTq25final[-(1:2),]
DPTq75final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
DPTq75final <- DPTq75final[-(1:2),]
DPTq90final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
DPTq90final <- DPTq90final[-(1:2),]
DPTminfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
DPTminfinal <- DPTminfinal[-(1:2),]
DPTmaxfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
DPTmaxfinal <- DPTmaxfinal[-(1:2),]
DPTHRESfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
DPTHRESfinal <- DPTHRESfinal[-(1:2),]


T2Mmeanfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
T2Mmeanfinal <- T2Mmeanfinal[-(1:2),]
T2Msdfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
T2Msdfinal <- T2Msdfinal[-(1:2),]
T2Mmedianfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
T2Mmedianfinal <- T2Mmedianfinal[-(1:2),]
T2Mq10final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
T2Mq10final <- T2Mq10final[-(1:2),]
T2Mq25final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
T2Mq25final <- T2Mq25final[-(1:2),]
T2Mq75final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
T2Mq75final <- T2Mq75final[-(1:2),]
T2Mq90final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
T2Mq90final <- T2Mq90final[-(1:2),]
T2Mminfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
T2Mminfinal <- T2Mminfinal[-(1:2),]
T2Mmaxfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
T2Mmaxfinal <- T2Mmaxfinal[-(1:2),]
T2MHRESfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
T2MHRESfinal <- T2MHRESfinal[-(1:2),]

EVAmeanfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
EVAmeanfinal <- EVAmeanfinal[-(1:2),]
EVAsdfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
EVAsdfinal <- EVAsdfinal[-(1:2),]
EVAmedianfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
EVAmedianfinal <- EVAmedianfinal[-(1:2),]
EVAq10final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
EVAq10final <- EVAq10final[-(1:2),]
EVAq25final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
EVAq25final <- EVAq25final[-(1:2),]
EVAq75final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
EVAq75final <- EVAq75final[-(1:2),]
EVAq90final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
EVAq90final <- EVAq90final[-(1:2),]
EVAminfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
EVAminfinal <- EVAminfinal[-(1:2),]
EVAmaxfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
EVAmaxfinal <- EVAmaxfinal[-(1:2),]
EVAHRESfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
EVAHRESfinal <- EVAHRESfinal[-(1:2),]


CCmeanfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CCmeanfinal <- CCmeanfinal[-(1:2),]
CCsdfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CCsdfinal <- CCsdfinal[-(1:2),]
CCmedianfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CCmedianfinal <- CCmedianfinal[-(1:2),]
CCq10final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CCq10final <- CCq10final[-(1:2),]
CCq25final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CCq25final <- CCq25final[-(1:2),]
CCq75final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CCq75final <- CCq75final[-(1:2),]
CCq90final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CCq90final <- CCq90final[-(1:2),]
CCminfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CCminfinal <- CCminfinal[-(1:2),]
CCmaxfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CCmaxfinal <- CCmaxfinal[-(1:2),]
CCHRESfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CCHRESfinal <- CCHRESfinal[-(1:2),]

CAPEmeanfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CAPEmeanfinal <- CAPEmeanfinal[-(1:2),]
CAPEsdfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CAPEsdfinal <- CAPEsdfinal[-(1:2),]
CAPEmedianfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CAPEmedianfinal <- CAPEmedianfinal[-(1:2),]
CAPEq10final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CAPEq10final <- CAPEq10final[-(1:2),]
CAPEq25final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CAPEq25final <- CAPEq25final[-(1:2),]
CAPEq75final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CAPEq75final <- CAPEq75final[-(1:2),]
CAPEq90final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CAPEq90final <- CAPEq90final[-(1:2),]
CAPEminfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CAPEminfinal <- CAPEminfinal[-(1:2),]
CAPEmaxfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CAPEmaxfinal <- CAPEmaxfinal[-(1:2),]
CAPEHRESfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CAPEHRESfinal <- CAPEHRESfinal[-(1:2),]

CAPESmeanfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CAPESmeanfinal <- CAPESmeanfinal[-(1:2),]
CAPESsdfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CAPESsdfinal <- CAPESsdfinal[-(1:2),]
CAPESmedianfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CAPESmedianfinal <- CAPESmedianfinal[-(1:2),]
CAPESq10final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CAPESq10final <- CAPESq10final[-(1:2),]
CAPESq25final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CAPESq25final <- CAPESq25final[-(1:2),]
CAPESq75final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CAPESq75final <- CAPESq75final[-(1:2),]
CAPESq90final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CAPESq90final <- CAPESq90final[-(1:2),]
CAPESminfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CAPESminfinal <- CAPESminfinal[-(1:2),]
CAPESmaxfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CAPESmaxfinal <- CAPESmaxfinal[-(1:2),]
CAPESHRESfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
CAPESHRESfinal <- CAPESHRESfinal[-(1:2),]


test = unique(RRmean$ECmatch.RmatchID)
for (i in test){
  #  subdata = Rdata[which(Rdata$Rtime == i+hours(24)),]   
  subdata = Rdata[which(Rdata$RmatchID == i+3),]
  
  #    rID = upscal_match[,1]
  #    dfID = upscal_match[,2]   
  #    subdata = rbind.data.frame(rID,dfID,subdata)   
  subdata = t(subdata)
  subdata = as.data.frame(subdata)
  colnames(subdata) = c(subdata[2,1])
  subRRmean = RRmean[which(RRmean$ECmatch.RmatchID == i),]
  subRRsd = RRsd[which(RRsd$ECmatch.RmatchID == i),]
  subRRmedian = RRmedian[which(RRmedian$ECmatch.RmatchID == i),]
  subRRq10 = RRq10[which(RRq10$ECmatch.RmatchID == i),]
  subRRq25 = RRq25[which(RRq25$ECmatch.RmatchID == i),]
  subRRq75 = RRq75[which(RRq75$ECmatch.RmatchID == i),]
  subRRq90 = RRq90[which(RRq90$ECmatch.RmatchID == i),]
  subRRp0.05 = RRp0.05[which(RRp0.05$ECmatch.RmatchID == i),]
  subRRp0.3 = RRp0.3[which(RRp0.3$ECmatch.RmatchID == i),]
  subRRp10 = RRp10[which(RRp10$ECmatch.RmatchID == i),]
  subRRp20 = RRp20[which(RRp20$ECmatch.RmatchID == i),]
  subRRmin = RRmin[which(RRmin$ECmatch.RmatchID == i),]
  subRRmax = RRmax[which(RRmax$ECmatch.RmatchID == i),]
  subRRHRES = RRHRES[which(RRHRES$ECmatch.RmatchID == i),]
  
  subCONRRmean = CONRRmean[which(CONRRmean$ECmatch.RmatchID == i),]
  subCONRRsd = CONRRsd[which(CONRRsd$ECmatch.RmatchID == i),]
  subCONRRmedian = CONRRmedian[which(CONRRmedian$ECmatch.RmatchID == i),]
  subCONRRq10 = CONRRq10[which(CONRRq10$ECmatch.RmatchID == i),]
  subCONRRq25 = CONRRq25[which(CONRRq25$ECmatch.RmatchID == i),]
  subCONRRq75 = CONRRq75[which(CONRRq75$ECmatch.RmatchID == i),]
  subCONRRq90 = CONRRq90[which(CONRRq90$ECmatch.RmatchID == i),]
  subCONRRp0.05 = CONRRp0.05[which(CONRRp0.05$ECmatch.RmatchID == i),]
  subCONRRp0.3 = CONRRp0.3[which(CONRRp0.3$ECmatch.RmatchID == i),]
  subCONRRp10 = CONRRp10[which(CONRRp10$ECmatch.RmatchID == i),]
  subCONRRp20 = CONRRp20[which(CONRRp20$ECmatch.RmatchID == i),]
  subCONRRmin = CONRRmin[which(CONRRmin$ECmatch.RmatchID == i),]
  subCONRRmax = CONRRmax[which(CONRRmax$ECmatch.RmatchID == i),]
  subCONRRHRES = CONRRHRES[which(CONRRHRES$ECmatch.RmatchID == i),]
  
  subWSPDmean = WSPDmean[which(WSPDmean$ECmatch.RmatchID == i),]
  subWSPDsd = WSPDsd[which(WSPDsd$ECmatch.RmatchID == i),]
  subWSPDmedian = WSPDmedian[which(WSPDmedian$ECmatch.RmatchID == i),]
  subWSPDq10 = WSPDq10[which(WSPDq10$ECmatch.RmatchID == i),]
  subWSPDq25 = WSPDq25[which(WSPDq25$ECmatch.RmatchID == i),]
  subWSPDq75 = WSPDq75[which(WSPDq75$ECmatch.RmatchID == i),]
  subWSPDq90 = WSPDq90[which(WSPDq90$ECmatch.RmatchID == i),]
  subWSPDmin = WSPDmin[which(WSPDmin$ECmatch.RmatchID == i),]
  subWSPDmax = WSPDmax[which(WSPDmax$ECmatch.RmatchID == i),]
  subWSPDHRES = WSPDHRES[which(WSPDHRES$ECmatch.RmatchID == i),]
  
  subDPTmean = DPTmean[which(DPTmean$ECmatch.RmatchID == i),]
  subDPTsd = DPTsd[which(DPTsd$ECmatch.RmatchID == i),]
  subDPTmedian = DPTmedian[which(DPTmedian$ECmatch.RmatchID == i),]
  subDPTq10 = DPTq10[which(DPTq10$ECmatch.RmatchID == i),]
  subDPTq25 = DPTq25[which(DPTq25$ECmatch.RmatchID == i),]
  subDPTq75 = DPTq75[which(DPTq75$ECmatch.RmatchID == i),]
  subDPTq90 = DPTq90[which(DPTq90$ECmatch.RmatchID == i),]
  subDPTmin = DPTmin[which(DPTmin$ECmatch.RmatchID == i),]
  subDPTmax = DPTmax[which(DPTmax$ECmatch.RmatchID == i),]
  subDPTHRES = DPTHRES[which(DPTHRES$ECmatch.RmatchID == i),]
  
  subT2Mmean = T2Mmean[which(T2Mmean$ECmatch.RmatchID == i),]
  subT2Msd = T2Msd[which(T2Msd$ECmatch.RmatchID == i),]
  subT2Mmedian = T2Mmedian[which(T2Mmedian$ECmatch.RmatchID== i),]
  subT2Mq10 = T2Mq10[which(T2Mq10$ECmatch.RmatchID== i),]
  subT2Mq25 = T2Mq25[which(T2Mq25$ECmatch.RmatchID== i),]
  subT2Mq75 = T2Mq75[which(T2Mq75$ECmatch.RmatchID== i),]
  subT2Mq90 = T2Mq90[which(T2Mq90$ECmatch.RmatchID== i),]
  subT2Mmin = T2Mmin[which(T2Mmin$ECmatch.RmatchID== i),]
  subT2Mmax = T2Mmax[which(T2Mmax$ECmatch.RmatchID== i),]
  subT2MHRES = T2MHRES[which(T2MHRES$ECmatch.RmatchID== i),]
  
  subEVAmean = EVAmean[which(EVAmean$ECmatch.RmatchID== i),]
  subEVAsd = EVAsd[which(EVAsd$ECmatch.RmatchID== i),]
  subEVAmedian = EVAmedian[which(EVAmedian$ECmatch.RmatchID== i),]
  subEVAq10 = EVAq10[which(EVAq10$ECmatch.RmatchID== i),]
  subEVAq25 = EVAq25[which(EVAq25$ECmatch.RmatchID== i),]
  subEVAq75 = EVAq75[which(EVAq75$ECmatch.RmatchID== i),]
  subEVAq90 = EVAq90[which(EVAq90$ECmatch.RmatchID== i),]
  subEVAmin = EVAmin[which(EVAmin$ECmatch.RmatchID== i),]
  subEVAmax = EVAmax[which(EVAmax$ECmatch.RmatchID== i),]
  subEVAHRES = EVAHRES[which(EVAHRES$ECmatch.RmatchID== i),]
  
  subCCmean = CCmean[which(CCmean$ECmatch.RmatchID== i),]
  subCCsd = CCsd[which(CCsd$ECmatch.RmatchID== i),]
  subCCmedian = CCmedian[which(CCmedian$ECmatch.RmatchID== i),]
  subCCq10 = CCq10[which(CCq10$ECmatch.RmatchID== i),]
  subCCq25 = CCq25[which(CCq25$ECmatch.RmatchID== i),]
  subCCq75 = CCq75[which(CCq75$ECmatch.RmatchID== i),]
  subCCq90 = CCq90[which(CCq90$ECmatch.RmatchID== i),]
  subCCmin = CCmin[which(CCmin$ECmatch.RmatchID== i),]
  subCCmax = CCmax[which(CCmax$ECmatch.RmatchID== i),]
  subCCHRES = CCHRES[which(CCHRES$ECmatch.RmatchID== i),]
  
  subCAPEmean = CAPEmean[which(CAPEmean$ECmatch.RmatchID== i),]
  subCAPEsd = CAPEsd[which(CAPEsd$ECmatch.RmatchID== i),]
  subCAPEmedian = CAPEmedian[which(CAPEmedian$ECmatch.RmatchID== i),]
  subCAPEq10 = CAPEq10[which(CAPEq10$ECmatch.RmatchID== i),]
  subCAPEq25 = CAPEq25[which(CAPEq25$ECmatch.RmatchID== i),]
  subCAPEq75 = CAPEq75[which(CAPEq75$ECmatch.RmatchID== i),]
  subCAPEq90 = CAPEq90[which(CAPEq90$ECmatch.RmatchID== i),]
  subCAPEmin = CAPEmin[which(CAPEmin$ECmatch.RmatchID== i),]
  subCAPEmax = CAPEmax[which(CAPEmax$ECmatch.RmatchID== i),]
  subCAPEHRES = CAPEHRES[which(CAPEHRES$ECmatch.RmatchID== i),]
  
  subCAPESmean = CAPESmean[which(CAPESmean$ECmatch.RmatchID== i),]
  subCAPESsd = CAPESsd[which(CAPESsd$ECmatch.RmatchID== i),]
  subCAPESmedian = CAPESmedian[which(CAPESmedian$ECmatch.RmatchID== i),]
  subCAPESq10 = CAPESq10[which(CAPESq10$ECmatch.RmatchID== i),]
  subCAPESq25 = CAPESq25[which(CAPESq25$ECmatch.RmatchID== i),]
  subCAPESq75 = CAPESq75[which(CAPESq75$ECmatch.RmatchID== i),]
  subCAPESq90 = CAPESq90[which(CAPESq90$ECmatch.RmatchID== i),]
  subCAPESmin = CAPESmin[which(CAPESmin$ECmatch.RmatchID== i),]
  subCAPESmax = CAPESmax[which(CAPESmax$ECmatch.RmatchID== i),]
  subCAPESHRES = CAPESHRES[which(CAPESHRES$ECmatch.RmatchID== i),]
  
  
  
  
  
  
  
  
  subpre = cbind.data.frame(subRRmean[,1],subRRsd[,1],subRRmedian[,1],subRRq10[,1],subRRq25[,1],subRRq75[,1],subRRq90[,1],subRRp0.05[,1],
                            subRRp0.3[,1],subRRp10[,1],subRRp20[,1],subRRmin[,1],subRRmax[,1],subRRHRES[,1],
                            subCONRRmean[,1],subCONRRsd[,1],subCONRRmedian[,1],subCONRRq10[,1],subCONRRq25[,1],subCONRRq75[,1],subCONRRq90[,1],subCONRRp0.05[,1],
                            subCONRRp0.3[,1],subCONRRp10[,1],subCONRRp20[,1],subCONRRmin[,1],subCONRRmax[,1],subCONRRHRES[,1],
                            subWSPDmean[,1],subWSPDsd[,1],subWSPDmedian[,1],subWSPDq10[,1],subWSPDq25[,1],subWSPDq75[,1],subWSPDq90[,1], 
                            subWSPDmin[,1],subWSPDmax[,1],subWSPDHRES[,1],
                            subDPTmean[,1],subDPTsd[,1],subDPTmedian[,1],subDPTq10[,1],subDPTq25[,1],subDPTq75[,1],subDPTq90[,1], 
                            subDPTmin[,1],subDPTmax[,1],subDPTHRES[,1],     
                            subT2Mmean[,1],subT2Msd[,1],subT2Mmedian[,1],subT2Mq10[,1],subT2Mq25[,1],subT2Mq75[,1],subT2Mq90[,1], 
                            subT2Mmin[,1],subT2Mmax[,1],subT2MHRES[,1],     
                            subEVAmean[,1],subEVAsd[,1],subEVAmedian[,1],subEVAq10[,1],subEVAq25[,1],subEVAq75[,1],subEVAq90[,1], 
                            subEVAmin[,1],subEVAmax[,1],subEVAHRES[,1],     
                            subCCmean[,1],subCCsd[,1],subCCmedian[,1],subCCq10[,1],subCCq25[,1],subCCq75[,1],subCCq90[,1], 
                            subCCmin[,1],subCCmax[,1],subCCHRES[,1],     
                            subCAPEmean[,1],subCAPEsd[,1],subCAPEmedian[,1],subCAPEq10[,1],subCAPEq25[,1],subCAPEq75[,1],subCAPEq90[,1], 
                            subCAPEmin[,1],subCAPEmax[,1],subCAPEHRES[,1],     
                            subCAPESmean[,1],subCAPESsd[,1],subCAPESmedian[,1],subCAPESq10[,1],subCAPESq25[,1],subCAPESq75[,1],subCAPESq90[,1], 
                            subCAPESmin[,1],subCAPESmax[,1],subCAPESHRES)   
  
  colnames(subpre)[1] = "subRRmean"
  colnames(subpre)[2] = "subRRsd"
  colnames(subpre)[3] = "subRRmedian"
  colnames(subpre)[4] = "subRRq10"
  colnames(subpre)[5] = "subRRq25"
  colnames(subpre)[6] = "subp0"
  colnames(subpre)[1:98] = c("subRRmean","subRRsd","subRRmedian","subRRq10","subRRq25","subRRq75","subRRq90","subRRp0.05",
                             "subRRp0.3","subRRp10","subRRp20","subRRmin","subRRmax","subRRHRES",
                             "subCONRRmean","subCONRRsd","subCONRRmedian","subCONRRq10","subCONRRq25","subCONRRq75","subCONRRq90","subCONRRp0.05",
                             "subCONRRp0.3","subCONRRp10","subCONRRp20","subCONRRmin","subCONRRmax","subCONRRHRES",
                             "subWSPDmean","subWSPDsd","subWSPDmedian","subWSPDq10","subWSPDq25","subWSPDq75","subWSPDq90",
                             "subWSPDmin","subWSPDmax","subWSPDHRES",
                             "subDPTmean","subDPTsd","subDPTmedian","subDPTq10","subDPTq25","subDPTq75","subDPTq90",
                             "subDPTmin","subDPTmax","subDPTHRES",
                             "subT2Mmean","subT2Msd","subT2Mmedian","subT2Mq10","subT2Mq25","subT2Mq75","subT2Mq90",
                             "subT2Mmin","subT2Mmax","subT2MHRES",
                             "subEVAmean","subEVAsd","subEVAmedian","subEVAq10","subEVAq25","subEVAq75","subEVAq90",
                             "subEVAmin","subEVAmax","subEVAHRES",
                             "subCCmean","subCCsd","subCCmedian","subCCq10","subCCq25","subCCq75","subCCq90",
                             "subCCmin","subCCmax","subCCHRES",
                             "subCAPEmean","subCAPEsd","subCAPEmedian","subCAPEq10","subCAPEq25","subCAPEq75","subCAPEq90",
                             "subCAPEmin","subCAPEmax","subCAPEHRES",
                             "subCAPESmean","subCAPESsd","subCAPESmedian","subCAPESq10","subCAPESq25","subCAPESq75","subCAPESq90",
                             "subCAPESmin","subCAPESmax","subCAPESHRES")
  
  
  subpre = subpre[,c(1:98,101,102)]
  subdata$dfID=as.numeric(upscal_match$dfID)
  subdata$rID=as.numeric(upscal_match$rID)
  match = left_join(subdata,subpre,by = c("dfID"="grid"))
  match = match[-(1:2),]
  ECRRmean = match[,4]
  ECRRsd = match[,5]
  ECRRmedian = match[,6]
  ECRRq10 = match[,7]
  ECRRq25 = match[,8]
  ECRRq75 = match[,9]
  ECRRq90 = match[,10]
  ECRRp0.05 = match[,11]
  ECRRp0.3 = match[,12]
  ECRRp10 = match[,13]
  ECRRp20 = match[,14]
  ECRRmin = match[,15]
  ECRRmax = match[,16]
  ECRRHRES = match[,17]
  ECCONRRmean = match[,18]
  ECCONRRsd = match[,19]
  ECCONRRmedian = match[,20]
  ECCONRRq10 = match[,21]
  ECCONRRq25 = match[,22]
  ECCONRRq75 = match[,23]
  ECCONRRq90 = match[,24]
  ECCONRRp0.05 = match[,25]
  ECCONRRp0.3 = match[,26]
  ECCONRRp10 = match[,27]
  ECCONRRp20 = match[,28]
  ECCONRRmin = match[,29]
  ECCONRRmax = match[,30]
  ECCONRRHRES = match[,31]
  ECWSPDmean = match[,32]
  ECWSPDsd = match[,33]
  ECWSPDmedian = match[,34]
  ECWSPDq10 = match[,35]
  ECWSPDq25 = match[,36]
  ECWSPDq75 = match[,37]
  ECWSPDq90 = match[,38]
  ECWSPDmin = match[,39]
  ECWSPDmax = match[,40]
  ECWSPDHRES = match[,41]
  ECDPTmean = match[,42]
  ECDPTsd = match[,43]
  ECDPTmedian = match[,44]
  ECDPTq10 = match[,45]
  ECDPTq25 = match[,46]
  ECDPTq75 = match[,47]
  ECDPTq90 = match[,48]
  ECDPTmin = match[,49]
  ECDPTmax = match[,50]
  ECDPTHRES = match[,51]
  ECT2Mmean = match[,52]
  ECT2Msd = match[,53]
  ECT2Mmedian = match[,51]
  ECT2Mq10 = match[,55]
  ECT2Mq25 = match[,56]
  ECT2Mq75 = match[,57]
  ECT2Mq90 = match[,58]
  ECT2Mmin = match[,59]
  ECT2Mmax = match[,60]
  ECT2MHRES = match[,61]
  ECEVAmean = match[,62]
  ECEVAsd = match[,63]
  ECEVAmedian = match[,64]
  ECEVAq10 = match[,65]
  ECEVAq25 = match[,66]
  ECEVAq75 = match[,67]
  ECEVAq90 = match[,68]
  ECEVAmin = match[,69]
  ECEVAmax = match[,70]
  ECEVAHRES = match[,71]
  ECCCmean = match[,72]
  ECCCsd = match[,73]
  ECCCmedian = match[,74]
  ECCCq10 = match[,75]
  ECCCq25 = match[,76]
  ECCCq75 = match[,77]
  ECCCq90 = match[,78]
  ECCCmin = match[,79]
  ECCCmax = match[,80]
  ECCCHRES = match[,81]
  ECCAPECmean = match[,82]
  ECCAPECsd = match[,83]
  ECCAPECmedian = match[,84]
  ECCAPECq10 = match[,85]
  ECCAPECq25 = match[,86]
  ECCAPECq75 = match[,87]
  ECCAPECq90 = match[,88]
  ECCAPECmin = match[,89]
  ECCAPECmax = match[,90]
  ECCAPECHRES = match[,91]
  ECCAPESCmean = match[,92]
  ECCAPESCsd = match[,93]
  ECCAPESCmedian = match[,94]
  ECCAPESCq10 = match[,95]
  ECCAPESCq25 = match[,96]
  ECCAPESCq75 = match[,97]
  ECCAPESCq90 = match[,98]
  ECCAPESCmin = match[,99]
  ECCAPESCmax = match[,100]
  ECCAPESCHRES = match[,101]
  R4match = match[,1]
  dftime = match[,102]
  ECID = match[,2]
  RID = match[,3]
  
  R4final = cbind.data.frame(R4final,R4match)
  EC4time <- cbind.data.frame(EC4time,dftime)
  EC4ID = cbind.data.frame(EC4ID,ECID)
  R4ID = cbind.data.frame(R4ID,RID)
  RRmeanfinal=cbind.data.frame(RRmeanfinal,ECRRmean)
  RRsdfinal=cbind.data.frame(RRsdfinal,ECRRsd)
  RRmedianfinal=cbind.data.frame(RRmedianfinal,ECRRmedian)
  RRq10final=cbind.data.frame(RRq10final,ECRRq10)
  RRq25final=cbind.data.frame(RRq25final,ECRRq25)
  RRq75final=cbind.data.frame(RRq75final,ECRRq75)
  RRq90final=cbind.data.frame(RRq90final,ECRRq90)
  RRP0.05final=cbind.data.frame(RRP0.05final,ECRRp0.05)
  RRP0.3final=cbind.data.frame(RRP0.3final,ECRRp0.3)
  RRP10final=cbind.data.frame(RRP10final,ECRRp10)
  RRP20final=cbind.data.frame(RRP20final,ECRRp20)
  RRminfinal=cbind.data.frame(RRminfinal,ECRRmin)
  RRmaxfinal=cbind.data.frame(RRmaxfinal,ECRRmax)
  RRHRESfinal=cbind.data.frame(RRHRESfinal,ECRRHRES)
  
  CONRRmeanfinal=cbind.data.frame(CONRRmeanfinal,ECCONRRmean)
  CONRRsdfinal=cbind.data.frame(CONRRsdfinal,ECCONRRsd)
  CONRRmedianfinal=cbind.data.frame(CONRRmedianfinal,ECCONRRmedian)
  CONRRq10final=cbind.data.frame(CONRRq10final,ECCONRRq10)
  CONRRq25final=cbind.data.frame(CONRRq25final,ECCONRRq25)
  CONRRq75final=cbind.data.frame(CONRRq75final,ECCONRRq75)
  CONRRq90final=cbind.data.frame(CONRRq90final,ECCONRRq90)
  CONRRP0.05final=cbind.data.frame(CONRRP0.05final,ECCONRRp0.05)
  CONRRP0.3final=cbind.data.frame(CONRRP0.3final,ECCONRRp0.3)
  CONRRP10final=cbind.data.frame(CONRRP10final,ECCONRRp10)
  CONRRP20final=cbind.data.frame(CONRRP20final,ECCONRRp20)
  CONRRminfinal=cbind.data.frame(CONRRminfinal,ECCONRRmin)
  CONRRmaxfinal=cbind.data.frame(CONRRmaxfinal,ECCONRRmax)
  CONRRHRESfinal=cbind.data.frame(CONRRHRESfinal,ECCONRRHRES)
  
  WSPDmeanfinal=cbind.data.frame(WSPDmeanfinal,ECWSPDmean)
  WSPDsdfinal=cbind.data.frame(WSPDsdfinal,ECWSPDsd)
  WSPDmedianfinal=cbind.data.frame(WSPDmedianfinal,ECWSPDmedian)
  WSPDq10final=cbind.data.frame(WSPDq10final,ECWSPDq10)
  WSPDq25final=cbind.data.frame(WSPDq25final,ECWSPDq25)
  WSPDq75final=cbind.data.frame(WSPDq75final,ECWSPDq75)
  WSPDq90final=cbind.data.frame(WSPDq90final,ECWSPDq90)
  WSPDminfinal=cbind.data.frame(WSPDminfinal,ECWSPDmin)
  WSPDmaxfinal=cbind.data.frame(WSPDmaxfinal,ECWSPDmax)
  WSPDHRESfinal=cbind.data.frame(WSPDHRESfinal,ECWSPDHRES)
  
  DPTmeanfinal=cbind.data.frame(DPTmeanfinal,ECDPTmean)
  DPTsdfinal=cbind.data.frame(DPTsdfinal,ECDPTsd)
  DPTmedianfinal=cbind.data.frame(DPTmedianfinal,ECDPTmedian)
  DPTq10final=cbind.data.frame(DPTq10final,ECDPTq10)
  DPTq25final=cbind.data.frame(DPTq25final,ECDPTq25)
  DPTq75final=cbind.data.frame(DPTq75final,ECDPTq75)
  DPTq90final=cbind.data.frame(DPTq90final,ECDPTq90)
  DPTminfinal=cbind.data.frame(DPTminfinal,ECDPTmin)
  DPTmaxfinal=cbind.data.frame(DPTmaxfinal,ECDPTmax)
  DPTHRESfinal=cbind.data.frame(DPTHRESfinal,ECDPTHRES)
  
  T2Mmeanfinal=cbind.data.frame(T2Mmeanfinal,ECT2Mmean)
  T2Msdfinal=cbind.data.frame(T2Msdfinal,ECT2Msd)
  T2Mmedianfinal=cbind.data.frame(T2Mmedianfinal,ECT2Mmedian)
  T2Mq10final=cbind.data.frame(T2Mq10final,ECT2Mq10)
  T2Mq25final=cbind.data.frame(T2Mq25final,ECT2Mq25)
  T2Mq75final=cbind.data.frame(T2Mq75final,ECT2Mq75)
  T2Mq90final=cbind.data.frame(T2Mq90final,ECT2Mq90)
  T2Mminfinal=cbind.data.frame(T2Mminfinal,ECT2Mmin)
  T2Mmaxfinal=cbind.data.frame(T2Mmaxfinal,ECT2Mmax)
  T2MHRESfinal=cbind.data.frame(T2MHRESfinal,ECT2MHRES)
  
  EVAmeanfinal=cbind.data.frame(EVAmeanfinal,ECEVAmean)
  EVAsdfinal=cbind.data.frame(EVAsdfinal,ECEVAsd)
  EVAmedianfinal=cbind.data.frame(EVAmedianfinal,ECEVAmedian)
  EVAq10final=cbind.data.frame(EVAq10final,ECEVAq10)
  EVAq25final=cbind.data.frame(EVAq25final,ECEVAq25)
  EVAq75final=cbind.data.frame(EVAq75final,ECEVAq75)
  EVAq90final=cbind.data.frame(EVAq90final,ECEVAq90)
  EVAminfinal=cbind.data.frame(EVAminfinal,ECEVAmin)
  EVAmaxfinal=cbind.data.frame(EVAmaxfinal,ECEVAmax)
  EVAHRESfinal=cbind.data.frame(EVAHRESfinal,ECEVAHRES)
  
  CCmeanfinal=cbind.data.frame(CCmeanfinal,ECCCmean)
  CCsdfinal=cbind.data.frame(CCsdfinal,ECCCsd)
  CCmedianfinal=cbind.data.frame(CCmedianfinal,ECCCmedian)
  CCq10final=cbind.data.frame(CCq10final,ECCCq10)
  CCq25final=cbind.data.frame(CCq25final,ECCCq25)
  CCq75final=cbind.data.frame(CCq75final,ECCCq75)
  CCq90final=cbind.data.frame(CCq90final,ECCCq90)
  CCminfinal=cbind.data.frame(CCminfinal,ECCCmin)
  CCmaxfinal=cbind.data.frame(CCmaxfinal,ECCCmax)
  CCHRESfinal=cbind.data.frame(CCHRESfinal,ECCCHRES)
  
  CAPEmeanfinal=cbind.data.frame(CAPEmeanfinal,ECCAPECmean)
  CAPEsdfinal=cbind.data.frame(CAPEsdfinal,ECCAPECsd)
  CAPEmedianfinal=cbind.data.frame(CAPEmedianfinal,ECCAPECmedian)
  CAPEq10final=cbind.data.frame(CAPEq10final,ECCAPECq10)
  CAPEq25final=cbind.data.frame(CAPEq25final,ECCAPECq25)
  CAPEq75final=cbind.data.frame(CAPEq75final,ECCAPECq75)
  CAPEq90final=cbind.data.frame(CAPEq90final,ECCAPECq90)
  CAPEminfinal=cbind.data.frame(CAPEminfinal,ECCAPECmin)
  CAPEmaxfinal=cbind.data.frame(CAPEmaxfinal,ECCAPECmax)
  CAPEHRESfinal=cbind.data.frame(CAPEHRESfinal,ECCAPECHRES)
  
  CAPESmeanfinal=cbind.data.frame(CAPESmeanfinal,ECCAPESCmean)
  CAPESsdfinal=cbind.data.frame(CAPESsdfinal,ECCAPESCsd)
  CAPESmedianfinal=cbind.data.frame(CAPESmedianfinal,ECCAPESCmedian)
  CAPESq10final=cbind.data.frame(CAPESq10final,ECCAPESCq10)
  CAPESq25final=cbind.data.frame(CAPESq25final,ECCAPESCq25)
  CAPESq75final=cbind.data.frame(CAPESq75final,ECCAPESCq75)
  CAPESq90final=cbind.data.frame(CAPESq90final,ECCAPESCq90)
  CAPESminfinal=cbind.data.frame(CAPESminfinal,ECCAPESCmin)
  CAPESmaxfinal=cbind.data.frame(CAPESmaxfinal,ECCAPESCmax)
  CAPESHRESfinal=cbind.data.frame(CAPESHRESfinal,ECCAPESCHRES)

}



RRmeanfinal1 = RRmeanfinal[,-1]
RRmeanfinal1 = as.matrix(RRmeanfinal1)
dim(RRmeanfinal1) <- c(599*602,1)
RRsdfinal1 = RRsdfinal[,-1]
RRsdfinal1 = as.matrix(RRsdfinal1)
dim(RRsdfinal1) <- c(599*602,1)
RRmedianfinal1 = RRmedianfinal[,-1]
RRmedianfinal1 = as.matrix(RRmedianfinal1)
dim(RRmedianfinal1) <- c(599*602,1)
RRq10final1 = RRq10final[,-1]
RRq10final1 = as.matrix(RRq10final1)
dim(RRq10final1) <- c(599*602,1)
RRq25final1 = RRq25final[,-1]
RRq25final1 = as.matrix(RRq25final1)
dim(RRq25final1) <- c(599*602,1)
RRq75final1 = RRq75final[,-1]
RRq75final1 = as.matrix(RRq75final1)
dim(RRq75final1) <- c(599*602,1)
RRq90final1 = RRq90final[,-1]
RRq90final1 = as.matrix(RRq90final1)
dim(RRq90final1) <- c(599*602,1)
RRP0.05final1 = RRP0.05final[,-1]
RRP0.05final1 = as.matrix(RRP0.05final1)
dim(RRP0.05final1) <- c(599*602,1)
RRP0.3final1 = RRP0.3final[,-1]
RRP0.3final1 = as.matrix(RRP0.3final1)
dim(RRP0.3final1) <- c(599*602,1)
RRP10final1 = RRP10final[,-1]
RRP10final1 = as.matrix(RRP10final1)
dim(RRP10final1) <- c(599*602,1)
RRP20final1 = RRP20final[,-1]
RRP20final1 = as.matrix(RRP20final1)
dim(RRP20final1) <- c(599*602,1)
RRminfinal1 = RRminfinal[,-1]
RRminfinal1 = as.matrix(RRminfinal1)
dim(RRminfinal1) <- c(599*602,1)
RRmaxfinal1 = RRmaxfinal[,-1]
RRmaxfinal1 = as.matrix(RRmaxfinal1)
dim(RRmaxfinal1) <- c(599*602,1)
RRHRESfinal1 = RRHRESfinal[,-1]
RRHRESfinal1 = as.matrix(RRHRESfinal1)
dim(RRHRESfinal1) <- c(599*602,1)

CONRRmeanfinal1 = CONRRmeanfinal[,-1]
CONRRmeanfinal1 = as.matrix(CONRRmeanfinal1)
dim(CONRRmeanfinal1) <- c(599*602,1)
CONRRsdfinal1 = CONRRsdfinal[,-1]
CONRRsdfinal1 = as.matrix(CONRRsdfinal1)
dim(CONRRsdfinal1) <- c(599*602,1)
CONRRmedianfinal1 = CONRRmedianfinal[,-1]
CONRRmedianfinal1 = as.matrix(CONRRmedianfinal1)
dim(CONRRmedianfinal1) <- c(599*602,1)
CONRRq10final1 = CONRRq10final[,-1]
CONRRq10final1 = as.matrix(CONRRq10final1)
dim(CONRRq10final1) <- c(599*602,1)
CONRRq25final1 = CONRRq25final[,-1]
CONRRq25final1 = as.matrix(CONRRq25final1)
dim(CONRRq25final1) <- c(599*602,1)
CONRRq75final1 = CONRRq75final[,-1]
CONRRq75final1 = as.matrix(CONRRq75final1)
dim(CONRRq75final1) <- c(599*602,1)
CONRRq90final1 = CONRRq90final[,-1]
CONRRq90final1 = as.matrix(CONRRq90final1)
dim(CONRRq90final1) <- c(599*602,1)
CONRRP0.05final1 = CONRRP0.05final[,-1]
CONRRP0.05final1 = as.matrix(CONRRP0.05final1)
dim(CONRRP0.05final1) <- c(599*602,1)
CONRRP0.3final1 = CONRRP0.3final[,-1]
CONRRP0.3final1 = as.matrix(CONRRP0.3final1)
dim(CONRRP0.3final1) <- c(599*602,1)
CONRRP10final1 = CONRRP10final[,-1]
CONRRP10final1 = as.matrix(CONRRP10final1)
dim(CONRRP10final1) <- c(599*602,1)
CONRRP20final1 = CONRRP20final[,-1]
CONRRP20final1 = as.matrix(CONRRP20final1)
dim(CONRRP20final1) <- c(599*602,1)
CONRRminfinal1 = CONRRminfinal[,-1]
CONRRminfinal1 = as.matrix(CONRRminfinal1)
dim(CONRRminfinal1) <- c(599*602,1)
CONRRmaxfinal1 = CONRRmaxfinal[,-1]
CONRRmaxfinal1 = as.matrix(CONRRmaxfinal1)
dim(CONRRmaxfinal1) <- c(599*602,1)
CONRRHRESfinal1 = CONRRHRESfinal[,-1]
CONRRHRESfinal1 = as.matrix(CONRRHRESfinal1)
dim(CONRRHRESfinal1) <- c(599*602,1)

WSPDmeanfinal1 = WSPDmeanfinal[,-1]
WSPDmeanfinal1 = as.matrix(WSPDmeanfinal1)
dim(WSPDmeanfinal1) <- c(599*602,1)
WSPDsdfinal1 = WSPDsdfinal[,-1]
WSPDsdfinal1 = as.matrix(WSPDsdfinal1)
dim(WSPDsdfinal1) <- c(599*602,1)
WSPDmedianfinal1 = WSPDmedianfinal[,-1]
WSPDmedianfinal1 = as.matrix(WSPDmedianfinal1)
dim(WSPDmedianfinal1) <- c(599*602,1)
WSPDq10final1 = WSPDq10final[,-1]
WSPDq10final1 = as.matrix(WSPDq10final1)
dim(WSPDq10final1) <- c(599*602,1)
WSPDq25final1 = WSPDq25final[,-1]
WSPDq25final1 = as.matrix(WSPDq25final1)
dim(WSPDq25final1) <- c(599*602,1)
WSPDq75final1 = WSPDq75final[,-1]
WSPDq75final1 = as.matrix(WSPDq75final1)
dim(WSPDq75final1) <- c(599*602,1)
WSPDq90final1 = WSPDq90final[,-1]
WSPDq90final1 = as.matrix(WSPDq90final1)
dim(WSPDq90final1) <- c(599*602,1)
WSPDminfinal1 = WSPDminfinal[,-1]
WSPDminfinal1 = as.matrix(WSPDminfinal1)
dim(WSPDminfinal1) <- c(599*602,1)
WSPDmaxfinal1 = WSPDmaxfinal[,-1]
WSPDmaxfinal1 = as.matrix(WSPDmaxfinal1)
dim(WSPDmaxfinal1) <- c(599*602,1)
WSPDHRESfinal1 = WSPDHRESfinal[,-1]
WSPDHRESfinal1 = as.matrix(WSPDHRESfinal1)
dim(WSPDHRESfinal1) <- c(599*602,1)  

DPTmeanfinal1 = DPTmeanfinal[,-1]
DPTmeanfinal1 = as.matrix(DPTmeanfinal1)
dim(DPTmeanfinal1) <- c(599*602,1)
DPTsdfinal1 = DPTsdfinal[,-1]
DPTsdfinal1 = as.matrix(DPTsdfinal1)
dim(DPTsdfinal1) <- c(599*602,1)
DPTmedianfinal1 = DPTmedianfinal[,-1]
DPTmedianfinal1 = as.matrix(DPTmedianfinal1)
dim(DPTmedianfinal1) <- c(599*602,1)
DPTq10final1 = DPTq10final[,-1]
DPTq10final1 = as.matrix(DPTq10final1)
dim(DPTq10final1) <- c(599*602,1)
DPTq25final1 = DPTq25final[,-1]
DPTq25final1 = as.matrix(DPTq25final1)
dim(DPTq25final1) <- c(599*602,1)
DPTq75final1 = DPTq75final[,-1]
DPTq75final1 = as.matrix(DPTq75final1)
dim(DPTq75final1) <- c(599*602,1)
DPTq90final1 = DPTq90final[,-1]
DPTq90final1 = as.matrix(DPTq90final1)
dim(DPTq90final1) <- c(599*602,1)
DPTminfinal1 = DPTminfinal[,-1]
DPTminfinal1 = as.matrix(DPTminfinal1)
dim(DPTminfinal1) <- c(599*602,1)
DPTmaxfinal1 = DPTmaxfinal[,-1]
DPTmaxfinal1 = as.matrix(DPTmaxfinal1)
dim(DPTmaxfinal1) <- c(599*602,1)
DPTHRESfinal1 = DPTHRESfinal[,-1]
DPTHRESfinal1 = as.matrix(DPTHRESfinal1)
dim(DPTHRESfinal1) <- c(599*602,1)  

T2Mmeanfinal1 = T2Mmeanfinal[,-1]
T2Mmeanfinal1 = as.matrix(T2Mmeanfinal1)
dim(T2Mmeanfinal1) <- c(599*602,1)
T2Msdfinal1 = T2Msdfinal[,-1]
T2Msdfinal1 = as.matrix(T2Msdfinal1)
dim(T2Msdfinal1) <- c(599*602,1)
T2Mmedianfinal1 = T2Mmedianfinal[,-1]
T2Mmedianfinal1 = as.matrix(T2Mmedianfinal1)
dim(T2Mmedianfinal1) <- c(599*602,1)
T2Mq10final1 = T2Mq10final[,-1]
T2Mq10final1 = as.matrix(T2Mq10final1)
dim(T2Mq10final1) <- c(599*602,1)
T2Mq25final1 = T2Mq25final[,-1]
T2Mq25final1 = as.matrix(T2Mq25final1)
dim(T2Mq25final1) <- c(599*602,1)
T2Mq75final1 = T2Mq75final[,-1]
T2Mq75final1 = as.matrix(T2Mq75final1)
dim(T2Mq75final1) <- c(599*602,1)
T2Mq90final1 = T2Mq90final[,-1]
T2Mq90final1 = as.matrix(T2Mq90final1)
dim(T2Mq90final1) <- c(599*602,1)
T2Mminfinal1 = T2Mminfinal[,-1]
T2Mminfinal1 = as.matrix(T2Mminfinal1)
dim(T2Mminfinal1) <- c(599*602,1)
T2Mmaxfinal1 = T2Mmaxfinal[,-1]
T2Mmaxfinal1 = as.matrix(T2Mmaxfinal1)
dim(T2Mmaxfinal1) <- c(599*602,1)
T2MHRESfinal1 = T2MHRESfinal[,-1]
T2MHRESfinal1 = as.matrix(T2MHRESfinal1)
dim(T2MHRESfinal1) <- c(599*602,1)  

EVAmeanfinal1 = EVAmeanfinal[,-1]
EVAmeanfinal1 = as.matrix(EVAmeanfinal1)
dim(EVAmeanfinal1) <- c(599*602,1)
EVAsdfinal1 = EVAsdfinal[,-1]
EVAsdfinal1 = as.matrix(EVAsdfinal1)
dim(EVAsdfinal1) <- c(599*602,1)
EVAmedianfinal1 = EVAmedianfinal[,-1]
EVAmedianfinal1 = as.matrix(EVAmedianfinal1)
dim(EVAmedianfinal1) <- c(599*602,1)
EVAq10final1 = EVAq10final[,-1]
EVAq10final1 = as.matrix(EVAq10final1)
dim(EVAq10final1) <- c(599*602,1)
EVAq25final1 = EVAq25final[,-1]
EVAq25final1 = as.matrix(EVAq25final1)
dim(EVAq25final1) <- c(599*602,1)
EVAq75final1 = EVAq75final[,-1]
EVAq75final1 = as.matrix(EVAq75final1)
dim(EVAq75final1) <- c(599*602,1)
EVAq90final1 = EVAq90final[,-1]
EVAq90final1 = as.matrix(EVAq90final1)
dim(EVAq90final1) <- c(599*602,1)
EVAminfinal1 = EVAminfinal[,-1]
EVAminfinal1 = as.matrix(EVAminfinal1)
dim(EVAminfinal1) <- c(599*602,1)
EVAmaxfinal1 = EVAmaxfinal[,-1]
EVAmaxfinal1 = as.matrix(EVAmaxfinal1)
dim(EVAmaxfinal1) <- c(599*602,1)
EVAHRESfinal1 = EVAHRESfinal[,-1]
EVAHRESfinal1 = as.matrix(EVAHRESfinal1)
dim(EVAHRESfinal1) <- c(599*602,1)  

CCmeanfinal1 = CCmeanfinal[,-1]
CCmeanfinal1 = as.matrix(CCmeanfinal1)
dim(CCmeanfinal1) <- c(599*602,1)
CCsdfinal1 = CCsdfinal[,-1]
CCsdfinal1 = as.matrix(CCsdfinal1)
dim(CCsdfinal1) <- c(599*602,1)
CCmedianfinal1 = CCmedianfinal[,-1]
CCmedianfinal1 = as.matrix(CCmedianfinal1)
dim(CCmedianfinal1) <- c(599*602,1)
CCq10final1 = CCq10final[,-1]
CCq10final1 = as.matrix(CCq10final1)
dim(CCq10final1) <- c(599*602,1)
CCq25final1 = CCq25final[,-1]
CCq25final1 = as.matrix(CCq25final1)
dim(CCq25final1) <- c(599*602,1)
CCq75final1 = CCq75final[,-1]
CCq75final1 = as.matrix(CCq75final1)
dim(CCq75final1) <- c(599*602,1)
CCq90final1 = CCq90final[,-1]
CCq90final1 = as.matrix(CCq90final1)
dim(CCq90final1) <- c(599*602,1)
CCminfinal1 = CCminfinal[,-1]
CCminfinal1 = as.matrix(CCminfinal1)
dim(CCminfinal1) <- c(599*602,1)
CCmaxfinal1 = CCmaxfinal[,-1]
CCmaxfinal1 = as.matrix(CCmaxfinal1)
dim(CCmaxfinal1) <- c(599*602,1)
CCHRESfinal1 = CCHRESfinal[,-1]
CCHRESfinal1 = as.matrix(CCHRESfinal1)
dim(CCHRESfinal1) <- c(599*602,1)  

CAPEmeanfinal1 = CAPEmeanfinal[,-1]
CAPEmeanfinal1 = as.matrix(CAPEmeanfinal1)
dim(CAPEmeanfinal1) <- c(599*602,1)
CAPEsdfinal1 = CAPEsdfinal[,-1]
CAPEsdfinal1 = as.matrix(CAPEsdfinal1)
dim(CAPEsdfinal1) <- c(599*602,1)
CAPEmedianfinal1 = CAPEmedianfinal[,-1]
CAPEmedianfinal1 = as.matrix(CAPEmedianfinal1)
dim(CAPEmedianfinal1) <- c(599*602,1)
CAPEq10final1 = CAPEq10final[,-1]
CAPEq10final1 = as.matrix(CAPEq10final1)
dim(CAPEq10final1) <- c(599*602,1)
CAPEq25final1 = CAPEq25final[,-1]
CAPEq25final1 = as.matrix(CAPEq25final1)
dim(CAPEq25final1) <- c(599*602,1)
CAPEq75final1 = CAPEq75final[,-1]
CAPEq75final1 = as.matrix(CAPEq75final1)
dim(CAPEq75final1) <- c(599*602,1)
CAPEq90final1 = CAPEq90final[,-1]
CAPEq90final1 = as.matrix(CAPEq90final1)
dim(CAPEq90final1) <- c(599*602,1)
CAPEminfinal1 = CAPEminfinal[,-1]
CAPEminfinal1 = as.matrix(CAPEminfinal1)
dim(CAPEminfinal1) <- c(599*602,1)
CAPEmaxfinal1 = CAPEmaxfinal[,-1]
CAPEmaxfinal1 = as.matrix(CAPEmaxfinal1)
dim(CAPEmaxfinal1) <- c(599*602,1)
CAPEHRESfinal1 = CAPEHRESfinal[,-1]
CAPEHRESfinal1 = as.matrix(CAPEHRESfinal1)
dim(CAPEHRESfinal1) <- c(599*602,1)  

CAPESmeanfinal1 = CAPESmeanfinal[,-1]
CAPESmeanfinal1 = as.matrix(CAPESmeanfinal1)
dim(CAPESmeanfinal1) <- c(599*602,1)
CAPESsdfinal1 = CAPESsdfinal[,-1]
CAPESsdfinal1 = as.matrix(CAPESsdfinal1)
dim(CAPESsdfinal1) <- c(599*602,1)
CAPESmedianfinal1 = CAPESmedianfinal[,-1]
CAPESmedianfinal1 = as.matrix(CAPESmedianfinal1)
dim(CAPESmedianfinal1) <- c(599*602,1)
CAPESq10final1 = CAPESq10final[,-1]
CAPESq10final1 = as.matrix(CAPESq10final1)
dim(CAPESq10final1) <- c(599*602,1)
CAPESq25final1 = CAPESq25final[,-1]
CAPESq25final1 = as.matrix(CAPESq25final1)
dim(CAPESq25final1) <- c(599*602,1)
CAPESq75final1 = CAPESq75final[,-1]
CAPESq75final1 = as.matrix(CAPESq75final1)
dim(CAPESq75final1) <- c(599*602,1)
CAPESq90final1 = CAPESq90final[,-1]
CAPESq90final1 = as.matrix(CAPESq90final1)
dim(CAPESq90final1) <- c(599*602,1)
CAPESminfinal1 = CAPESminfinal[,-1]
CAPESminfinal1 = as.matrix(CAPESminfinal1)
dim(CAPESminfinal1) <- c(599*602,1)
CAPESmaxfinal1 = CAPESmaxfinal[,-1]
CAPESmaxfinal1 = as.matrix(CAPESmaxfinal1)
dim(CAPESmaxfinal1) <- c(599*602,1)
CAPESHRESfinal1 = CAPESHRESfinal[,-1]
CAPESHRESfinal1 = as.matrix(CAPESHRESfinal1)
dim(CAPESHRESfinal1) <- c(599*602,1)  



R4final1 = R4final[,-1]
R4final1 = as.matrix(R4final1)
dim(R4final1) <- c(599*602,1)
EC4time1 = EC4time[,-1]
EC4time1 = as.matrix(EC4time1)
dim(EC4time1) <- c(599*602,1)
EC4ID1 = EC4ID[,-1]
EC4ID1 = as.matrix(EC4ID1)
dim(EC4ID1) <- c(599*602,1)
R4ID1 = R4ID[,-1]
R4ID1 = as.matrix(R4ID1)
dim(R4ID1) <- c(599*602,1)


final = cbind.data.frame(RRmeanfinal1,RRsdfinal1, RRmedianfinal1, RRq10final1,RRq25final1,RRq75final1,RRq90final1, 
                         RRP0.05final1,RRP0.3final1,RRP10final1,RRP20final1,RRminfinal1,RRmaxfinal1,RRHRESfinal1,
                         CONRRmeanfinal1,CONRRsdfinal1,CONRRmedianfinal1,CONRRq10final1,CONRRq25final1,
                         CONRRq75final1,CONRRq90final1,CONRRP0.05final1,CONRRP0.3final1,CONRRP10final1,
                         CONRRP20final1,CONRRminfinal1,CONRRmaxfinal1,CONRRHRESfinal1,
                         WSPDmeanfinal1,WSPDsdfinal1,WSPDmedianfinal1,WSPDq10final1,WSPDq25final1,WSPDq75final1,
                         WSPDq90final1,WSPDminfinal1,WSPDmaxfinal1,WSPDHRESfinal1,
                         DPTmeanfinal1,DPTsdfinal1,DPTmedianfinal1,DPTq10final1,DPTq25final1,DPTq75final1,
                         DPTq90final1,DPTminfinal1,DPTmaxfinal1,DPTHRESfinal1,
                         T2Mmeanfinal1,T2Msdfinal1,T2Mmedianfinal1,T2Mq10final1,T2Mq25final1,T2Mq75final1,
                         T2Mq90final1,T2Mminfinal1,T2Mmaxfinal1,T2MHRESfinal1,
                         EVAmeanfinal1,EVAsdfinal1,EVAmedianfinal1,EVAq10final1,EVAq25final1,EVAq75final1,
                         EVAq90final1,EVAminfinal1,EVAmaxfinal1,EVAHRESfinal1,
                         CCmeanfinal1,CCsdfinal1,CCmedianfinal1,CCq10final1,CCq25final1,CCq75final1,
                         CCq90final1,CCminfinal1,CCmaxfinal1,CCHRESfinal1,
                         CAPEmeanfinal1,CAPEsdfinal1,CAPEmedianfinal1,CAPEq10final1,CAPEq25final1,CAPEq75final1,
                         CAPEq90final1,CAPEminfinal1,CAPEmaxfinal1,CAPEHRESfinal1,
                         CAPESmeanfinal1,CAPESsdfinal1,CAPESmedianfinal1,CAPESq10final1,CAPESq25final1,CAPESq75final1,
                         CAPESq90final1,CAPESminfinal1,CAPESmaxfinal1,CAPESHRESfinal1,
                         R4final1,EC4time1,EC4ID1,R4ID1)

final[which(final$RRmeanfinal1 <0 ),1]<- 0
#setwd("E:/output")
#write.table(final, "24hour_forecast_alldata.txt", sep = "\t", col.names = T, row.names = F)
final = read.table("E:/output/24hour_forecast_alldata.txt", sep = "\t", header = T)

final = na.omit(final)
final$R4final1 = as.numeric(final$R4final1)


#########Make stratified random sample############################################################

library(gamlss)
library(dplyr)
set.seed(1)
final1 <- final %>%
  group_by(EC4time1) %>%
  sample_n(300)

final2 = setdiff(final, final1)
#######divide the training and validation by year#################################################
final2018 = final[1:164126,]
final2020 = setdiff(final,final2018)
########Model with ZAGA#########################################################################
mod0<-gamlss(R4final1~1, data=final2020, family=ZAGA)

mod1 = stepGAICAll.A(mod0, scope=list(lower=~1,upper=~RRmeanfinal1+RRsdfinal1+ RRmedianfinal1+ RRq10final1+RRq25final1+RRq75final1+RRq90final1+ 
                                        RRP0.05final1+RRP0.3final1+RRP10final1+RRP20final1+RRminfinal1+RRmaxfinal1+RRHRESfinal1+
                                        CONRRmeanfinal1+CONRRsdfinal1+CONRRmedianfinal1+CONRRq10final1+CONRRq25final1+
                                        CONRRq75final1+CONRRq90final1+CONRRP0.05final1+CONRRP0.3final1+CONRRP10final1+
                                        CONRRP20final1+CONRRminfinal1+CONRRmaxfinal1+CONRRHRESfinal1+
                                        WSPDmeanfinal1+WSPDsdfinal1+WSPDmedianfinal1+WSPDq10final1+WSPDq25final1+WSPDq75final1+
                                        WSPDq90final1+WSPDminfinal1+WSPDmaxfinal1+WSPDHRESfinal1+
                                        DPTmeanfinal1+DPTsdfinal1+DPTmedianfinal1+DPTq10final1+DPTq25final1+DPTq75final1+
                                        DPTq90final1+DPTminfinal1+DPTmaxfinal1+DPTHRESfinal1+
                                        T2Mmeanfinal1+T2Msdfinal1+T2Mmedianfinal1+T2Mq10final1+T2Mq25final1+T2Mq75final1+
                                        T2Mq90final1+T2Mminfinal1+T2Mmaxfinal1+T2MHRESfinal1+
                                        EVAmeanfinal1+EVAsdfinal1+EVAmedianfinal1+EVAq10final1+EVAq25final1+EVAq75final1+
                                        EVAq90final1+EVAminfinal1+EVAmaxfinal1+EVAHRESfinal1+
                                        CCmeanfinal1+CCsdfinal1+CCmedianfinal1+CCq10final1+CCq25final1+CCq75final1+
                                        CCq90final1+CCminfinal1+CCmaxfinal1+CCHRESfinal1+
                                        CAPEmeanfinal1+CAPEsdfinal1+CAPEmedianfinal1+CAPEq10final1+CAPEq25final1+CAPEq75final1+
                                        CAPEq90final1+CAPEminfinal1+CAPEmaxfinal1+CAPEHRESfinal1+
                                        CAPESmeanfinal1+CAPESsdfinal1+CAPESmedianfinal1+CAPESq10final1+CAPESq25final1+CAPESq75final1+
                                        CAPESq90final1+CAPESminfinal1+CAPESmaxfinal1+CAPESHRESfinal1), 
                     sigma.scope=list(lower=~1,upper=~RRmeanfinal1+RRsdfinal1+ RRmedianfinal1+ RRq10final1+RRq25final1+RRq75final1+RRq90final1+ 
                                        RRP0.05final1+RRP0.3final1+RRP10final1+RRP20final1+RRminfinal1+RRmaxfinal1+RRHRESfinal1+
                                        CONRRmeanfinal1+CONRRsdfinal1+CONRRmedianfinal1+CONRRq10final1+CONRRq25final1+
                                        CONRRq75final1+CONRRq90final1+CONRRP0.05final1+CONRRP0.3final1+CONRRP10final1+
                                        CONRRP20final1+CONRRminfinal1+CONRRmaxfinal1+CONRRHRESfinal1+
                                        WSPDmeanfinal1+WSPDsdfinal1+WSPDmedianfinal1+WSPDq10final1+WSPDq25final1+WSPDq75final1+
                                        WSPDq90final1+WSPDminfinal1+WSPDmaxfinal1+WSPDHRESfinal1+
                                        DPTmeanfinal1+DPTsdfinal1+DPTmedianfinal1+DPTq10final1+DPTq25final1+DPTq75final1+
                                        DPTq90final1+DPTminfinal1+DPTmaxfinal1+DPTHRESfinal1+
                                        T2Mmeanfinal1+T2Msdfinal1+T2Mmedianfinal1+T2Mq10final1+T2Mq25final1+T2Mq75final1+
                                        T2Mq90final1+T2Mminfinal1+T2Mmaxfinal1+T2MHRESfinal1+
                                        EVAmeanfinal1+EVAsdfinal1+EVAmedianfinal1+EVAq10final1+EVAq25final1+EVAq75final1+
                                        EVAq90final1+EVAminfinal1+EVAmaxfinal1+EVAHRESfinal1+
                                        CCmeanfinal1+CCsdfinal1+CCmedianfinal1+CCq10final1+CCq25final1+CCq75final1+
                                        CCq90final1+CCminfinal1+CCmaxfinal1+CCHRESfinal1+
                                        CAPEmeanfinal1+CAPEsdfinal1+CAPEmedianfinal1+CAPEq10final1+CAPEq25final1+CAPEq75final1+
                                        CAPEq90final1+CAPEminfinal1+CAPEmaxfinal1+CAPEHRESfinal1+
                                        CAPESmeanfinal1+CAPESsdfinal1+CAPESmedianfinal1+CAPESq10final1+CAPESq25final1+CAPESq75final1+
                                        CAPESq90final1+CAPESminfinal1+CAPESmaxfinal1+CAPESHRESfinal1), 
                     nu.scope=list(lower=~1,upper=~RRmeanfinal1+RRsdfinal1+ RRmedianfinal1+ RRq10final1+RRq25final1+RRq75final1+RRq90final1+ 
                                     RRP0.05final1+RRP0.3final1+RRP10final1+RRP20final1+RRminfinal1+RRmaxfinal1+RRHRESfinal1+
                                     CONRRmeanfinal1+CONRRsdfinal1+CONRRmedianfinal1+CONRRq10final1+CONRRq25final1+
                                     CONRRq75final1+CONRRq90final1+CONRRP0.05final1+CONRRP0.3final1+CONRRP10final1+
                                     CONRRP20final1+CONRRminfinal1+CONRRmaxfinal1+CONRRHRESfinal1+
                                     WSPDmeanfinal1+WSPDsdfinal1+WSPDmedianfinal1+WSPDq10final1+WSPDq25final1+WSPDq75final1+
                                     WSPDq90final1+WSPDminfinal1+WSPDmaxfinal1+WSPDHRESfinal1+
                                     DPTmeanfinal1+DPTsdfinal1+DPTmedianfinal1+DPTq10final1+DPTq25final1+DPTq75final1+
                                     DPTq90final1+DPTminfinal1+DPTmaxfinal1+DPTHRESfinal1+
                                     T2Mmeanfinal1+T2Msdfinal1+T2Mmedianfinal1+T2Mq10final1+T2Mq25final1+T2Mq75final1+
                                     T2Mq90final1+T2Mminfinal1+T2Mmaxfinal1+T2MHRESfinal1+
                                     EVAmeanfinal1+EVAsdfinal1+EVAmedianfinal1+EVAq10final1+EVAq25final1+EVAq75final1+
                                     EVAq90final1+EVAminfinal1+EVAmaxfinal1+EVAHRESfinal1+
                                     CCmeanfinal1+CCsdfinal1+CCmedianfinal1+CCq10final1+CCq25final1+CCq75final1+
                                     CCq90final1+CCminfinal1+CCmaxfinal1+CCHRESfinal1+
                                     CAPEmeanfinal1+CAPEsdfinal1+CAPEmedianfinal1+CAPEq10final1+CAPEq25final1+CAPEq75final1+
                                     CAPEq90final1+CAPEminfinal1+CAPEmaxfinal1+CAPEHRESfinal1+
                                     CAPESmeanfinal1+CAPESsdfinal1+CAPESmedianfinal1+CAPESq10final1+CAPESq25final1+CAPESq75final1+
                                     CAPESq90final1+CAPESminfinal1+CAPESmaxfinal1+CAPESHRESfinal1),
                     #direction = "forward",
                     steps = 4) 

######Make the predict for the new data##########################################################################
prd  <- predictAll(mod1,newdata=final2018)
prob <- 1:51 / 52
new_data = as.data.frame(matrix(NA,nrow = 51, ncol = dim(final2018)[1]))
for (i in 1:dim(final2018)[1]){
  new_data[,i]=sapply(prob, function(u) qZAGA(u, mu=prd$mu[i], sigma=prd$sigma[i], nu=prd$nu[i]))
}
new_data = as.data.frame(t(new_data))
######Calculate CRPS for each group of the data
new_data$observation = final2018$R4final1
library(scoringRules)
crps = as.data.frame(matrix(NA,nrow = dim(final2018)[1],ncol = 1))

for (i in 1:dim(final2018)[1]){
  sample = as.numeric(new_data[i,1:51])
  crps[i,]=crps_sample(y = new_data$observation[i],dat = sample)
}
#setwd("E:/output")
write.table(crps, "crps_24hours_4steps_new.txt", sep = "\t", col.names = T, row.names = F)
crps = read.table("E:/output/crps_24hours_3steps_new.txt", sep = "\t", header = T)
