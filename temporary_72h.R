#For calculate the CDF
CDF = as.data.frame(matrix(NA,nrow = 15, ncol = dim(final2018)[1]))
for (i in 1:dim(final2018)[1]){
  CDF[1,i]=sapply(0.05, function(u) pZAGA(u, mu=prd$mu[i], sigma=prd$sigma[i], nu=prd$nu[i]))
  CDF[2,i]=sapply(0.1, function(u) pZAGA(u, mu=prd$mu[i], sigma=prd$sigma[i], nu=prd$nu[i]))
  CDF[3,i]=sapply(0.2, function(u) pZAGA(u, mu=prd$mu[i], sigma=prd$sigma[i], nu=prd$nu[i]))
  CDF[4,i]=sapply(0.5, function(u) pZAGA(u, mu=prd$mu[i], sigma=prd$sigma[i], nu=prd$nu[i]))
  CDF[5,i]=sapply(0.8, function(u) pZAGA(u, mu=prd$mu[i], sigma=prd$sigma[i], nu=prd$nu[i]))
  CDF[6,i]=sapply(1, function(u) pZAGA(u, mu=prd$mu[i], sigma=prd$sigma[i], nu=prd$nu[i]))
  CDF[7,i]=sapply(2, function(u) pZAGA(u, mu=prd$mu[i], sigma=prd$sigma[i], nu=prd$nu[i]))
  CDF[8,i]=sapply(3, function(u) pZAGA(u, mu=prd$mu[i], sigma=prd$sigma[i], nu=prd$nu[i]))
  CDF[9,i]=sapply(5, function(u) pZAGA(u, mu=prd$mu[i], sigma=prd$sigma[i], nu=prd$nu[i]))
  CDF[10,i]=sapply(7.5, function(u) pZAGA(u, mu=prd$mu[i], sigma=prd$sigma[i], nu=prd$nu[i]))
  CDF[11,i]=sapply(10, function(u) pZAGA(u, mu=prd$mu[i], sigma=prd$sigma[i], nu=prd$nu[i]))
  CDF[12,i]=sapply(12.5, function(u) pZAGA(u, mu=prd$mu[i], sigma=prd$sigma[i], nu=prd$nu[i]))
  CDF[13,i]=sapply(15, function(u) pZAGA(u, mu=prd$mu[i], sigma=prd$sigma[i], nu=prd$nu[i]))
  CDF[14,i]=sapply(17.5, function(u) pZAGA(u, mu=prd$mu[i], sigma=prd$sigma[i], nu=prd$nu[i]))
  CDF[15,i]=sapply(20, function(u) pZAGA(u, mu=prd$mu[i], sigma=prd$sigma[i], nu=prd$nu[i]))
}
CDF = as.data.frame(t(CDF))
write.table(CDF, "NCDF_72hour_3steps.txt", sep = "\t", col.names = T, row.names = F)
crps3step = mean(crps3$crps)
rm(CDF,mod0,mod1,clim,new_data,prd,crps3)





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
saveRDS(object = mod1, file = "72h_forecast_3rd_year_4step.rds")




prd  <- predictAll(mod1,newdata=final2018)
prob <- 1:51 / 52
new_data = as.data.frame(matrix(NA,nrow = 51, ncol = dim(final2018)[1]))
for (i in 1:dim(final2018)[1]){
  new_data[,i]=sapply(prob, function(u) qZAGA(u, mu=prd$mu[i], sigma=prd$sigma[i], nu=prd$nu[i]))
}
new_data = as.data.frame(t(new_data))
#write.table(new_data, "new_data_72hours_3steps.txt", sep = "\t", col.names = T, row.names = F)

new_data$observation = final2018$R4final1
library(scoringRules)
crps = as.data.frame(matrix(NA,nrow = dim(final2018)[1],ncol = 1))

for (i in 1:dim(final2018)[1]){
  sample = as.numeric(new_data[i,1:51])
  crps[i,]=crps_sample(y = new_data$observation[i],dat = sample)
}
crps4step = mean(crps[,1])
rm(mod0,mod1,new_data,prd,crps)



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
                     steps = 5) 
saveRDS(object = mod1, file = "72h_forecast_3rd_year_5step.rds")




prd  <- predictAll(mod1,newdata=final2018)
prob <- 1:51 / 52
new_data = as.data.frame(matrix(NA,nrow = 51, ncol = dim(final2018)[1]))
for (i in 1:dim(final2018)[1]){
  new_data[,i]=sapply(prob, function(u) qZAGA(u, mu=prd$mu[i], sigma=prd$sigma[i], nu=prd$nu[i]))
}
new_data = as.data.frame(t(new_data))
#write.table(new_data, "new_data_72hours_3steps.txt", sep = "\t", col.names = T, row.names = F)

new_data$observation = final2018$R4final1
library(scoringRules)
crps = as.data.frame(matrix(NA,nrow = dim(final2018)[1],ncol = 1))

for (i in 1:dim(final2018)[1]){
  sample = as.numeric(new_data[i,1:51])
  crps[i,]=crps_sample(y = new_data$observation[i],dat = sample)
}
crps5step = mean(crps[,1])
rm(mod0,mod1,new_data,prd,crps)