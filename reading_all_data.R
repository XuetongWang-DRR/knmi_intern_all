
#Only select the fourth column of forecast data (for forecasting after 19-24 hours)
RRmean = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_mean_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_mean_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_mean_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_mean_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_mean_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_mean_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_mean_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_mean_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_mean_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_mean_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_mean_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_mean_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_mean_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



RRsd = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_sd_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_sd_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_sd_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_sd_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_sd_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_sd_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_sd_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_sd_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_sd_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_sd_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_sd_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_sd_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_sd_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

RRmedian = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_median_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_median_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_median_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_median_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_median_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_median_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_median_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_median_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_median_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_median_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_median_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_median_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_median_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


RRq10 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_q10_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_q10_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q10_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q10_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q10_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q10_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q10_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q10_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q10_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q10_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q10_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q10_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q10_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

RRq25 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_q25_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_q25_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q25_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q25_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q25_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q25_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q25_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q25_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q25_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q25_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q25_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q25_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q25_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

RRq75 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_q75_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_q75_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q75_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q75_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q75_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q75_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q75_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q75_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q75_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q75_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q75_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q75_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q75_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


RRq90 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_q90_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_q90_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q90_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q90_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q90_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q90_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q90_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q90_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q90_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q90_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q90_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q90_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q90_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

RRp0.05 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_p0.05_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_p0.05_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p0.05_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p0.05_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p0.05_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p0.05_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_p0.05_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_p0.05_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_p0.05_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_p0.05_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_p0.05_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_p0.05_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_p0.05_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


RRp0.3 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_p0.3_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_p0.3_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p0.3_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p0.3_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p0.3_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p0.3_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_p0.3_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_p0.3_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_p0.3_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_p0.3_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_p0.3_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_p0.3_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_p0.3_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



RRp10 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_p10_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_p10_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p10_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p10_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p10_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p10_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_p10_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_p10_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_p10_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_p10_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_p10_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_p10_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_p10_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



RRp20 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_p20_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_p20_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p20_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p20_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p20_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p20_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_p20_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_p20_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_p20_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_p20_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_p20_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_p20_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_p20_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])




RRmin = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_min_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_min_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_min_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_min_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_min_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_min_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_min_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_min_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_min_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_min_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_min_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_min_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_min_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



RRmax = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_max_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_max_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_max_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_max_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_max_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_max_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_max_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_max_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_max_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_max_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_max_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_max_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_max_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



RRHRES = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_HRES_2018_11.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_HRES_2018_12.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_HRES_2019_01.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_HRES_2019_02.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_HRES_2019_03.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_HRES_2019_04.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_HRES_2020_10.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_HRES_2020_11.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_HRES_2020_12.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_HRES_2021_01.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_HRES_2021_02.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_HRES_2021_03.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_HRES_2021_04.txt",sep = "", header = T)[-1,c(4,41,42)])


WSPDmean = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/WSPD_6h_mean_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2018/WSPD_6h_mean_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_mean_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_mean_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_mean_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_mean_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_mean_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_mean_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_mean_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_mean_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_mean_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_mean_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_mean_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])





WSPDsd = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/WSPD_6h_sd_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/WSPD_6h_sd_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_sd_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_sd_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_sd_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_sd_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_sd_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_sd_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_sd_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_sd_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_sd_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_sd_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_sd_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


WSPDmedian = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/WSPD_6h_median_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2018/WSPD_6h_median_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_median_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_median_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_median_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_median_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_median_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_median_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_median_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_median_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_median_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_median_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_median_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])





WSPDq10 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/WSPD_6h_q10_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2018/WSPD_6h_q10_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_q10_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_q10_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_q10_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_q10_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_q10_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_q10_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_q10_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_q10_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_q10_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_q10_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_q10_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

WSPDq25 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/WSPD_6h_q25_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2018/WSPD_6h_q25_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_q25_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_q25_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_q25_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_q25_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_q25_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_q25_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_q25_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_q25_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_q25_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_q25_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_q25_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

WSPDq75 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/WSPD_6h_q75_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2018/WSPD_6h_q75_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_q75_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_q75_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_q75_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_q75_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_q75_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_q75_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_q75_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_q75_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_q75_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_q75_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_q75_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


WSPDq90 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/WSPD_6h_q90_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2018/WSPD_6h_q90_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_q90_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_q90_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_q90_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_q90_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_q90_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_q90_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_q90_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_q90_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_q90_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_q90_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_q90_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


WSPDmin = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/WSPD_6h_min_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2018/WSPD_6h_min_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_min_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_min_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_min_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_min_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_min_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_min_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_min_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_min_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_min_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_min_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_min_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



WSPDmax = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/WSPD_6h_max_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2018/WSPD_6h_max_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_max_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_max_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_max_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_max_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_max_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_max_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_max_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_max_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_max_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_max_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_max_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



WSPDHRES = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/WSPD_6h_HRES_2018_11.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2018/WSPD_6h_HRES_2018_12.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_HRES_2019_01.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_HRES_2019_02.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_HRES_2019_03.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/WSPD_6h_HRES_2019_04.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_HRES_2020_10.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_HRES_2020_11.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/WSPD_6h_HRES_2020_12.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_HRES_2021_01.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_HRES_2021_02.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_HRES_2021_03.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/WSPD_6h_HRES_2021_04.txt",sep = "", header = T)[-1,c(4,41,42)])






T2Mmean = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/T2M_6h_mean_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2018/T2M_6h_mean_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_mean_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_mean_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_mean_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_mean_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_mean_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_mean_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_mean_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_mean_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_mean_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_mean_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_mean_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

T2Msd = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/T2M_6h_sd_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2018/T2M_6h_sd_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_sd_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_sd_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_sd_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_sd_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_sd_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_sd_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_sd_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_sd_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_sd_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_sd_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_sd_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


T2Mmedian = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/T2M_6h_median_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2018/T2M_6h_median_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_median_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_median_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_median_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_median_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_median_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_median_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_median_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_median_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_median_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_median_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_median_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])





T2Mq10 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/T2M_6h_q10_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/T2M_6h_q10_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_q10_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_q10_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_q10_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_q10_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_q10_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_q10_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_q10_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_q10_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_q10_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_q10_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_q10_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

T2Mq25 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/T2M_6h_q25_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/T2M_6h_q25_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_q25_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_q25_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_q25_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_q25_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_q25_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_q25_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_q25_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_q25_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_q25_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_q25_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_q25_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

T2Mq75 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/T2M_6h_q75_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/T2M_6h_q75_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_q75_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_q75_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_q75_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_q75_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_q75_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_q75_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_q75_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_q75_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_q75_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_q75_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_q75_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


T2Mq90 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/T2M_6h_q90_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/T2M_6h_q90_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_q90_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_q90_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_q90_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_q90_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_q90_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_q90_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_q90_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_q90_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_q90_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_q90_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_q90_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


T2Mmin = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/T2M_6h_min_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/T2M_6h_min_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_min_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_min_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_min_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_min_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_min_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_min_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_min_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_min_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_min_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_min_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_min_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



T2Mmax = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/T2M_6h_max_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/T2M_6h_max_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_max_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_max_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_max_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_max_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_max_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_max_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_max_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_max_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_max_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_max_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_max_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



T2MHRES = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/T2M_6h_HRES_2018_11.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2018/T2M_6h_HRES_2018_12.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_HRES_2019_01.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_HRES_2019_02.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_HRES_2019_03.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/T2M_6h_HRES_2019_04.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_HRES_2020_10.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_HRES_2020_11.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/T2M_6h_HRES_2020_12.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_HRES_2021_01.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_HRES_2021_02.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_HRES_2021_03.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/T2M_6h_HRES_2021_04.txt",sep = "", header = T)[-1,c(4,41,42)])



EVAmean = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/EVA_6h_mean_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2018/EVA_6h_mean_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_mean_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_mean_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_mean_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_mean_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_mean_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_mean_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_mean_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_mean_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_mean_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_mean_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_mean_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

EVAsd = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/EVA_6h_sd_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2018/EVA_6h_sd_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_sd_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_sd_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_sd_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_sd_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_sd_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_sd_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_sd_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_sd_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_sd_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_sd_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_sd_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


EVAmedian = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/EVA_6h_median_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2018/EVA_6h_median_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_median_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_median_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_median_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_median_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_median_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_median_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_median_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_median_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_median_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_median_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_median_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])





EVAq10 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/EVA_6h_q10_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/EVA_6h_q10_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_q10_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_q10_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_q10_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_q10_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_q10_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_q10_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_q10_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_q10_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_q10_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_q10_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_q10_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

EVAq25 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/EVA_6h_q25_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/EVA_6h_q25_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_q25_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_q25_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_q25_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_q25_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_q25_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_q25_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_q25_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_q25_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_q25_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_q25_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_q25_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

EVAq75 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/EVA_6h_q75_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/EVA_6h_q75_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_q75_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_q75_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_q75_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_q75_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_q75_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_q75_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_q75_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_q75_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_q75_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_q75_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_q75_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


EVAq90 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/EVA_6h_q90_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/EVA_6h_q90_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_q90_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_q90_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_q90_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_q90_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_q90_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_q90_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_q90_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_q90_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_q90_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_q90_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_q90_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


EVAmin = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/EVA_6h_min_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/EVA_6h_min_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_min_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_min_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_min_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_min_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_min_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_min_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_min_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_min_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_min_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_min_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_min_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



EVAmax = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/EVA_6h_max_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/EVA_6h_max_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_max_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_max_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_max_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_max_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_max_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_max_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_max_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_max_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_max_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_max_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_max_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



EVAHRES = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/EVA_6h_HRES_2018_11.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2018/EVA_6h_HRES_2018_12.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_HRES_2019_01.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_HRES_2019_02.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_HRES_2019_03.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/EVA_6h_HRES_2019_04.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_HRES_2020_10.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_HRES_2020_11.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/EVA_6h_HRES_2020_12.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_HRES_2021_01.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_HRES_2021_02.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_HRES_2021_03.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/EVA_6h_HRES_2021_04.txt",sep = "", header = T)[-1,c(4,41,42)])




DPTmean = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/DPT_6h_mean_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2018/DPT_6h_mean_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_mean_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_mean_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_mean_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_mean_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_mean_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_mean_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_mean_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_mean_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_mean_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_mean_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_mean_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

DPTsd = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/DPT_6h_sd_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2018/DPT_6h_sd_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_sd_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_sd_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_sd_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_sd_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_sd_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_sd_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_sd_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_sd_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_sd_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_sd_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_sd_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


DPTmedian = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/DPT_6h_median_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2018/DPT_6h_median_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_median_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_median_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_median_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_median_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_median_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_median_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_median_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_median_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_median_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_median_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_median_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])





DPTq10 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/DPT_6h_q10_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/DPT_6h_q10_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_q10_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_q10_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_q10_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_q10_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_q10_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_q10_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_q10_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_q10_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_q10_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_q10_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_q10_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

DPTq25 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/DPT_6h_q25_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/DPT_6h_q25_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_q25_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_q25_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_q25_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_q25_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_q25_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_q25_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_q25_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_q25_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_q25_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_q25_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_q25_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

DPTq75 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/DPT_6h_q75_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/DPT_6h_q75_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_q75_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_q75_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_q75_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_q75_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_q75_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_q75_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_q75_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_q75_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_q75_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_q75_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_q75_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


DPTq90 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/DPT_6h_q90_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/DPT_6h_q90_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_q90_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_q90_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_q90_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_q90_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_q90_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_q90_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_q90_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_q90_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_q90_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_q90_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_q90_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


DPTmin = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/DPT_6h_min_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/DPT_6h_min_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_min_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_min_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_min_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_min_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_min_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_min_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_min_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_min_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_min_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_min_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_min_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



DPTmax = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/DPT_6h_max_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/DPT_6h_max_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_max_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_max_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_max_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_max_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_max_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_max_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_max_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_max_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_max_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_max_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_max_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



DPTHRES = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/DPT_6h_HRES_2018_11.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2018/DPT_6h_HRES_2018_12.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_HRES_2019_01.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_HRES_2019_02.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_HRES_2019_03.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/DPT_6h_HRES_2019_04.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_HRES_2020_10.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_HRES_2020_11.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/DPT_6h_HRES_2020_12.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_HRES_2021_01.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_HRES_2021_02.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_HRES_2021_03.txt",sep = "", header = T)[-1,c(4,41,42)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/DPT_6h_HRES_2021_04.txt",sep = "", header = T)[-1,c(4,41,42)])




CONRRmean = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_mean_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_mean_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_mean_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_mean_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_mean_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_mean_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_mean_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_mean_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_mean_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_mean_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_mean_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_mean_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_mean_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

CONRRsd = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_sd_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_sd_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_sd_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_sd_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_sd_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_sd_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_sd_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_sd_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_sd_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_sd_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_sd_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_sd_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_sd_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


CONRRmedian = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_median_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_median_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_median_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_median_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_median_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_median_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_median_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_median_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_median_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_median_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_median_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_median_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_median_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])





CONRRq10 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_q10_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_q10_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_q10_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_q10_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_q10_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_q10_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_q10_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_q10_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_q10_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_q10_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_q10_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_q10_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_q10_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

CONRRq25 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_q25_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_q25_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_q25_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_q25_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_q25_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_q25_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_q25_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_q25_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_q25_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_q25_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_q25_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_q25_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_q25_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

CONRRq75 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_q75_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_q75_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_q75_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_q75_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_q75_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_q75_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_q75_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_q75_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_q75_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_q75_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_q75_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_q75_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_q75_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


CONRRq90 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_q90_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_q90_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_q90_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_q90_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_q90_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_q90_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_q90_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_q90_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_q90_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_q90_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_q90_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_q90_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_q90_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

CONRRp0.05 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_p0.05_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_p0.05_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_p0.05_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_p0.05_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_p0.05_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_p0.05_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_p0.05_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_p0.05_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_p0.05_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_p0.05_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_p0.05_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_p0.05_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_p0.05_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


CONRRp0.3 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_p0.3_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_p0.3_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_p0.3_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_p0.3_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_p0.3_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_p0.3_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_p0.3_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_p0.3_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_p0.3_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_p0.3_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_p0.3_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_p0.3_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_p0.3_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



CONRRp10 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_p10_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_p10_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_p10_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_p10_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_p10_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_p10_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_p10_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_p10_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_p10_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_p10_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_p10_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_p10_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_p10_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



CONRRp20 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_p20_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_p20_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_p20_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_p20_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_p20_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_p20_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_p20_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_p20_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_p20_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_p20_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_p20_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_p20_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_p20_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])




CONRRmin = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_min_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_min_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_min_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_min_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_min_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_min_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_min_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_min_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_min_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_min_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_min_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_min_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_min_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



CONRRmax = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_max_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_max_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_max_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_max_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_max_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_max_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_max_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_max_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_max_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_max_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_max_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_max_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_max_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



CONRRHRES = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_HRES_2018_11.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2018/CONRR_6h_HRES_2018_12.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_HRES_2019_01.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_HRES_2019_02.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_HRES_2019_03.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/CONRR_6h_HRES_2019_04.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_HRES_2020_10.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_HRES_2020_11.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2020/CONRR_6h_HRES_2020_12.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_HRES_2021_01.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_HRES_2021_02.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_HRES_2021_03.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/CONRR_6h_HRES_2021_04.txt",sep = "", header = T)[-1,c(4,41,42)])




CCmean = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CC_6h_mean_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/CC_6h_mean_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_mean_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_mean_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_mean_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_mean_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_mean_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_mean_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_mean_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_mean_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_mean_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_mean_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_mean_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

CCsd = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CC_6h_sd_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2018/CC_6h_sd_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_sd_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_sd_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_sd_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_sd_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_sd_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_sd_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_sd_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_sd_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_sd_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_sd_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                        read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_sd_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


CCmedian = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CC_6h_median_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2018/CC_6h_median_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_median_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_median_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_median_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_median_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_median_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_median_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_median_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_median_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_median_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_median_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_median_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])





CCq10 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CC_6h_q10_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2018/CC_6h_q10_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_q10_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_q10_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_q10_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_q10_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_q10_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_q10_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_q10_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_q10_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_q10_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_q10_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_q10_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

CCq25 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CC_6h_q25_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2018/CC_6h_q25_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_q25_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_q25_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_q25_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_q25_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_q25_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_q25_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_q25_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_q25_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_q25_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_q25_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_q25_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

CCq75 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CC_6h_q75_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2018/CC_6h_q75_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_q75_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_q75_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_q75_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_q75_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_q75_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_q75_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_q75_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_q75_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_q75_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_q75_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_q75_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


CCq90 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CC_6h_q90_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2018/CC_6h_q90_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_q90_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_q90_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_q90_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_q90_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_q90_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_q90_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_q90_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_q90_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_q90_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_q90_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_q90_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


CCmin = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CC_6h_min_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2018/CC_6h_min_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_min_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_min_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_min_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_min_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_min_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_min_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_min_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_min_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_min_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_min_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_min_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



CCmax = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CC_6h_max_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2018/CC_6h_max_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_max_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_max_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_max_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_max_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_max_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_max_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_max_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_max_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_max_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_max_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                         read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_max_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



CCHRES = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CC_6h_HRES_2018_11.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/CC_6h_HRES_2018_12.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_HRES_2019_01.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_HRES_2019_02.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_HRES_2019_03.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/CC_6h_HRES_2019_04.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_HRES_2020_10.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_HRES_2020_11.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/CC_6h_HRES_2020_12.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_HRES_2021_01.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_HRES_2021_02.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_HRES_2021_03.txt",sep = "", header = T)[-1,c(4,41,42)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/CC_6h_HRES_2021_04.txt",sep = "", header = T)[-1,c(4,41,42)])





CAPEmean = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CAPE_6h_mean_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2018/CAPE_6h_mean_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_mean_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_mean_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_mean_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_mean_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_mean_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_mean_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_mean_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_mean_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_mean_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_mean_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_mean_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

CAPEsd = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CAPE_6h_sd_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2018/CAPE_6h_sd_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_sd_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_sd_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_sd_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_sd_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_sd_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_sd_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_sd_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_sd_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_sd_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_sd_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                          read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_sd_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


CAPEmedian = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CAPE_6h_median_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2018/CAPE_6h_median_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_median_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_median_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_median_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_median_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_median_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_median_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_median_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_median_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_median_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_median_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                              read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_median_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])





CAPEq10 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CAPE_6h_q10_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2018/CAPE_6h_q10_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_q10_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_q10_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_q10_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_q10_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_q10_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_q10_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_q10_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_q10_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_q10_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_q10_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_q10_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

CAPEq25 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CAPE_6h_q25_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2018/CAPE_6h_q25_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_q25_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_q25_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_q25_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_q25_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_q25_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_q25_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_q25_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_q25_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_q25_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_q25_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_q25_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

CAPEq75 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CAPE_6h_q75_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2018/CAPE_6h_q75_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_q75_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_q75_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_q75_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_q75_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_q75_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_q75_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_q75_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_q75_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_q75_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_q75_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_q75_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


CAPEq90 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CAPE_6h_q90_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2018/CAPE_6h_q90_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_q90_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_q90_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_q90_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_q90_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_q90_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_q90_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_q90_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_q90_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_q90_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_q90_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_q90_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


CAPEmin = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CAPE_6h_min_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2018/CAPE_6h_min_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_min_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_min_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_min_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_min_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_min_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_min_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_min_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_min_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_min_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_min_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_min_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



CAPEmax = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CAPE_6h_max_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2018/CAPE_6h_max_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_max_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_max_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_max_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_max_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_max_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_max_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_max_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_max_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_max_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_max_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_max_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



CAPEHRES = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CAPE_6h_HRES_2018_11.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2018/CAPE_6h_HRES_2018_12.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_HRES_2019_01.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_HRES_2019_02.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_HRES_2019_03.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPE_6h_HRES_2019_04.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_HRES_2020_10.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_HRES_2020_11.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPE_6h_HRES_2020_12.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_HRES_2021_01.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_HRES_2021_02.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_HRES_2021_03.txt",sep = "", header = T)[-1,c(4,41,42)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPE_6h_HRES_2021_04.txt",sep = "", header = T)[-1,c(4,41,42)])




CAPESmean = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CAPES_6h_mean_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2018/CAPES_6h_mean_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_mean_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_mean_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_mean_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_mean_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_mean_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_mean_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_mean_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_mean_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_mean_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_mean_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_mean_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

CAPESsd = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CAPES_6h_sd_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2018/CAPES_6h_sd_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_sd_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_sd_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_sd_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_sd_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_sd_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_sd_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_sd_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_sd_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_sd_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_sd_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                           read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_sd_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


CAPESmedian = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CAPES_6h_median_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2018/CAPES_6h_median_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_median_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_median_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_median_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_median_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_median_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_median_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_median_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_median_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_median_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_median_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                               read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_median_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])





CAPESq10 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CAPES_6h_q10_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2018/CAPES_6h_q10_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_q10_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_q10_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_q10_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_q10_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_q10_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_q10_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_q10_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_q10_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_q10_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_q10_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_q10_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

CAPESq25 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CAPES_6h_q25_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2018/CAPES_6h_q25_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_q25_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_q25_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_q25_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_q25_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_q25_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_q25_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_q25_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_q25_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_q25_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_q25_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_q25_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])

CAPESq75 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CAPES_6h_q75_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2018/CAPES_6h_q75_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_q75_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_q75_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_q75_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_q75_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_q75_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_q75_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_q75_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_q75_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_q75_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_q75_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_q75_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


CAPESq90 = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CAPES_6h_q90_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2018/CAPES_6h_q90_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_q90_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_q90_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_q90_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_q90_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_q90_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_q90_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_q90_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_q90_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_q90_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_q90_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_q90_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])


CAPESmin = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CAPES_6h_min_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2018/CAPES_6h_min_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_min_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_min_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_min_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_min_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_min_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_min_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_min_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_min_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_min_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_min_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_min_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



CAPESmax = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CAPES_6h_max_2018_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2018/CAPES_6h_max_2018_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_max_2019_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_max_2019_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_max_2019_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_max_2019_04.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_max_2020_10.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_max_2020_11.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_max_2020_12.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_max_2021_01.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_max_2021_02.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_max_2021_03.txt",sep = "", header = T)[-1,c(4,61,62)],
                            read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_max_2021_04.txt",sep = "", header = T)[-1,c(4,61,62)])



CAPESHRES = rbind.data.frame(read.table("E:/ECMWFdata/pre-process/Y2018/CAPES_6h_HRES_2018_11.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2018/CAPES_6h_HRES_2018_12.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_HRES_2019_01.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_HRES_2019_02.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_HRES_2019_03.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2019/CAPES_6h_HRES_2019_04.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_HRES_2020_10.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_HRES_2020_11.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2020/CAPES_6h_HRES_2020_12.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_HRES_2021_01.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_HRES_2021_02.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_HRES_2021_03.txt",sep = "", header = T)[-1,c(4,41,42)],
                             read.table("E:/ECMWFdata/pre-process/Y2021/CAPES_6h_HRES_2021_04.txt",sep = "", header = T)[-1,c(4,41,42)])


Rdata = cbind.data.frame(read.table("E:/output/radar_output/2018_11_9km_avg.txt",sep = "", header = T),
                         read.table("E:/output/radar_output/2018_12_9km_avg.txt",sep = "", header = T),
                         read.table("E:/output/radar_output/2019_01_9km_avg.txt",sep = "", header = T),
                         read.table("E:/output/radar_output/2019_02_9km_avg.txt",sep = "", header = T),
                         read.table("E:/output/radar_output/2019_03_9km_avg.txt",sep = "", header = T),
                         read.table("E:/output/radar_output/2019_04_9km_avg.txt",sep = "", header = T),
                         read.table("E:/output/radar_output/2019_10_9km_avg.txt",sep = "", header = T),
                         read.table("E:/output/radar_output/2019_11_9km_avg.txt",sep = "", header = T),
                         read.table("E:/output/radar_output/2019_12_9km_avg.txt",sep = "", header = T),
                         read.table("E:/output/radar_output/2020_01_9km_avg.txt",sep = "", header = T),
                         read.table("E:/output/radar_output/2020_02_9km_avg.txt",sep = "", header = T),
                         read.table("E:/output/radar_output/2020_03_9km_avg.txt",sep = "", header = T),
                         read.table("E:/output/radar_output/2020_04_9km_avg.txt",sep = "", header = T),
                         read.table("E:/output/radar_output/2020_10_9km_avg.txt",sep = "", header = T),
                         read.table("E:/output/radar_output/2020_11_9km_avg.txt",sep = "", header = T),
                         read.table("E:/output/radar_output/2020_12_9km_avg.txt",sep = "", header = T),
                         read.table("E:/output/radar_output/2021_01_9km_avg.txt",sep = "", header = T),
                         read.table("E:/output/radar_output/2021_02_9km_avg.txt",sep = "", header = T),
                         read.table("E:/output/radar_output/2021_03_9km_avg.txt",sep = "", header = T),
                         read.table("E:/output/radar_output/2021_04_9km_avg.txt",sep = "", header = T))



