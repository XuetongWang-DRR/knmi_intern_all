####parameter tunning#######################################################
crps_24h_raw = read.table("E:/output/validating_model/24hours/crps_24h_raw_forecast.txt", sep = "\t", header = T)
crps_24h_zaga = read.table("E:/output/validating_model/24hours/crps_24h_3steps.txt", sep = "\t", header = T)
crps_24h_qrf = read.table("E:/output/validating_model/24hours/crps_24h_ALLsample_5_3_500.txt", sep = "\t", header = T)



crps_48h_raw = read.table("E:/output/validating_model/48hours/crps_48h_raw_forecast.txt", sep = "\t", header = T)
crps_48h_zaga = read.table("E:/output/validating_model/48hours/crps_48h_3steps.txt", sep = "\t", header = T)
crps_48h_qrf = read.table("E:/output/validating_model/48hours/crps_48h_ALLsample_5_3_500.txt", sep = "\t", header = T)

crps_24h_zaga_2steps = read.table("E:/output/other_data/24hours/crps_24hours_2steps_new.txt", sep = "\t", header = T)
colnames(crps_24h_zaga_2steps) = c("crps")
crps_24h_zaga_2steps$crps_clim = crps_24h_zaga$crps_clim
crps_24h_zaga_2steps$crpss=(mean(crps_24h_zaga_2steps$crps)-mean(crps_24h_zaga_2steps$crps_clim))/(-mean(crps_24h_zaga_2steps$crps_clim))


crps_24h_zaga_4steps = read.table("E:/output/other_data/24hours/crps_24hours_4steps_new.txt", sep = "\t", header = T)
colnames(crps_24h_zaga_4steps) = c("crps")
crps_24h_zaga_4steps$crps_clim = crps_24h_zaga$crps_clim
crps_24h_zaga_4steps$crpss=(mean(crps_24h_zaga_4steps$crps)-mean(crps_24h_zaga_4steps$crps_clim))/(-mean(crps_24h_zaga_4steps$crps_clim))



crps_24h_zaga_5steps = read.table("E:/output/other_data/24hours/crps_24hours_5steps_new.txt", sep = "\t", header = T)
colnames(crps_24h_zaga_5steps) = c("crps")
crps_24h_zaga_5steps$crps_clim = crps_24h_zaga$crps_clim
crps_24h_zaga_5steps$crpss=(mean(crps_24h_zaga_5steps$crps)-mean(crps_24h_zaga_5steps$crps_clim))/(-mean(crps_24h_zaga_5steps$crps_clim))



crps_24h_zaga_6steps = read.table("E:/output/other_data/24hours/crps_24hours_6steps_new.txt", sep = "\t", header = T)
colnames(crps_24h_zaga_6steps) = c("crps")
crps_24h_zaga_6steps$crps_clim = crps_24h_zaga$crps_clim
crps_24h_zaga_6steps$crpss=(mean(crps_24h_zaga_6steps$crps)-mean(crps_24h_zaga_6steps$crps_clim))/(-mean(crps_24h_zaga_6steps$crps_clim))



crps_24h_zaga_7steps = read.table("E:/output/other_data/24hours/crps_24hours_7steps_new.txt", sep = "\t", header = T)
colnames(crps_24h_zaga_7steps) = c("crps")
crps_24h_zaga_7steps$crps_clim = crps_24h_zaga$crps_clim
crps_24h_zaga_7steps$crpss=(mean(crps_24h_zaga_7steps$crps)-mean(crps_24h_zaga_7steps$crps_clim))/(-mean(crps_24h_zaga_7steps$crps_clim))




crps_24h_zaga_8steps = read.table("E:/output/other_data/24hours/crps_24hours_8steps_new.txt", sep = "\t", header = T)
colnames(crps_24h_zaga_8steps) = c("crps")
crps_24h_zaga_8steps$crps_clim = crps_24h_zaga$crps_clim
crps_24h_zaga_8steps$crpss=(mean(crps_24h_zaga_8steps$crps)-mean(crps_24h_zaga_8steps$crps_clim))/(-mean(crps_24h_zaga_8steps$crps_clim))

crpss_24h_zaga_2steps = mean(crps_24h_zaga_2steps$crpss)
crpss_24h_zaga_3steps = mean(crps_24h_zaga$crpss)
crpss_24h_zaga_4steps = mean(crps_24h_zaga_4steps$crpss)
crpss_24h_zaga_5steps = mean(crps_24h_zaga_5steps$crpss)
crpss_24h_zaga_6steps = mean(crps_24h_zaga_6steps$crpss)
crpss_24h_zaga_7steps = mean(crps_24h_zaga_7steps$crpss)
crpss_24h_zaga_8steps = mean(crps_24h_zaga_8steps$crpss)

crpss_24h_raw = mean(crps_24h_raw$crpss)


crpss_24h_50sample_5_3_500 = (0.2637851-0.487895)/-0.487895
crpss_24h_100sample_5_3_500 = (0.2644968-0.487895)/-0.487895
crpss_24h_200sample_5_3_500 = (0.2640345-0.487895)/-0.487895
crpss_24h_ALLsample_5_3_500 = (0.2632237-0.487895)/-0.487895

crpss_24h_100sample_1_3_500 = (0.2643077-0.487895)/-0.487895
crpss_24h_100sample_2_3_500 = (0.2640682-0.487895)/-0.487895
crpss_24h_100sample_3_3_500 = (0.2640758-0.487895)/-0.487895
crpss_24h_100sample_4_3_500 = (0.2648345-0.487895)/-0.487895
crpss_24h_100sample_10_3_500 = (0.2650477-0.487895)/-0.487895

crpss_24h_100sample_5_4_500 = (0.2633883-0.487895)/-0.487895
crpss_24h_100sample_5_5_500 = (0.2625866-0.487895)/-0.487895
crpss_24h_100sample_5_6_500 = (0.2631154-0.487895)/-0.487895
crpss_24h_100sample_5_7_500 = (0.2625162-0.487895)/-0.487895
crpss_24h_100sample_5_8_500 = (0.262413-0.487895)/-0.487895
crpss_24h_100sample_5_9_500 = (0.262723-0.487895)/-0.487895
crpss_24h_100sample_5_10_500 = (0.2623146-0.487895)/-0.487895

crpss_24h_100sample_5_3_750 = (0.2645484-0.487895)/-0.487895
crpss_24h_100sample_5_3_1000 = (0.2641434-0.487895)/-0.487895
library(ggplot2)

#crpss针对zaga的对比
x <- c("2steps","3steps","4steps","5steps","6steps","7steps","8steps")

y <- c(crpss_24h_zaga_2steps,crpss_24h_zaga_3steps,crpss_24h_zaga_4steps,
       crpss_24h_zaga_5steps,crpss_24h_zaga_6steps,crpss_24h_zaga_7steps,
       crpss_24h_zaga_8steps)
df <- data.frame(x = x, y = y)
ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point(size = 5)+xlab("steps")+ylab("crpss")+
  ggtitle("crpss comparison for ZAGA with 24h forecast lead time") 


x = rep(c("24h"),times = 26)
y = c(crpss_24h_zaga_2steps,crpss_24h_zaga_3steps,crpss_24h_zaga_4steps,
      crpss_24h_zaga_5steps,crpss_24h_zaga_6steps,crpss_24h_zaga_7steps,
      crpss_24h_zaga_8steps,crpss_24h_50sample_5_3_500,crpss_24h_100sample_5_3_500,
      crpss_24h_200sample_5_3_500,crpss_24h_ALLsample_5_3_500,crpss_24h_100sample_1_3_500,
      crpss_24h_100sample_2_3_500,crpss_24h_100sample_3_3_500,crpss_24h_100sample_4_3_500,
      crpss_24h_100sample_10_3_500,crpss_24h_100sample_5_4_500,
      crpss_24h_100sample_5_5_500,crpss_24h_100sample_5_6_500,crpss_24h_100sample_5_7_500,
      crpss_24h_100sample_5_8_500,crpss_24h_100sample_5_9_500,crpss_24h_100sample_5_10_500,
      crpss_24h_100sample_5_3_750,crpss_24h_100sample_5_3_1000,crpss_24h_raw)

z = c("zaga","zaga","zaga","zaga","zaga","zaga","zaga","qrf","qrf","qrf","qrf","qrf","qrf","qrf","qrf"
      ,"qrf","qrf","qrf","qrf","qrf","qrf","qrf","qrf","qrf","qrf","raw")
df <- data.frame(x = x, y = y, z = z)
df$z = factor(df$z)

ggplot(data = df, mapping = aes(x = x, y = y, colour = z)) + geom_point(size = 3)+xlab("forecast time")+ylab("crpss")+
  ggtitle("crpss comparison") 





x <- c("Z1","Z2","Z3","Z4","Z5","Z6","Z7","RAW","Q1","Q2","Q3","Q4",
       "Q5","Q6","Q7","Q8","Q9","Q10","Q11","Q12","Q13","Q14",
       "Q15","Q16","Q17","Q18")
y <-  c(crpss_24h_zaga_2steps,crpss_24h_zaga_3steps,crpss_24h_zaga_4steps,
        crpss_24h_zaga_5steps,crpss_24h_zaga_6steps,crpss_24h_zaga_7steps,
        crpss_24h_zaga_8steps,crpss_24h_raw,crpss_24h_50sample_5_3_500,crpss_24h_100sample_5_3_500,
        crpss_24h_200sample_5_3_500,crpss_24h_ALLsample_5_3_500,crpss_24h_100sample_1_3_500,
        crpss_24h_100sample_2_3_500,crpss_24h_100sample_3_3_500,crpss_24h_100sample_4_3_500,
        crpss_24h_100sample_10_3_500,crpss_24h_100sample_5_4_500,
        crpss_24h_100sample_5_5_500,crpss_24h_100sample_5_6_500,crpss_24h_100sample_5_7_500,
        crpss_24h_100sample_5_8_500,crpss_24h_100sample_5_9_500,crpss_24h_100sample_5_10_500,
        crpss_24h_100sample_5_3_750,crpss_24h_100sample_5_3_1000)
z = c("zaga","zaga","zaga","zaga","zaga","zaga","zaga","raw","qrf","qrf","qrf","qrf","qrf","qrf","qrf","qrf"
      ,"qrf","qrf","qrf","qrf","qrf","qrf","qrf","qrf","qrf","qrf")
df <- data.frame(x = x, y = y, z = z)
#将数值型变量转换为因子型变量
df$z <- factor(df$z)
df$x = factor(df$x,levels=unique(df$x))
#分组变量赋值给颜色属性
ggplot(data = df, mapping = aes(x = x, y = y, colour = z)) + geom_point(size = 3)+
  scale_color_manual(values = c("blue","black","red"),name = "Type")+
  xlab("Model")+ylab("crpss")+ggtitle("crpss comparison")+
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))






