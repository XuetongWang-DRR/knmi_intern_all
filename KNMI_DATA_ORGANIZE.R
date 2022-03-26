####match two datasets------------------------------------------------------------------------------------------
library(ncdf4)
library(raster)
library(RANN)
library(gamlss.dist)  #ZAGA
data = read.delim("E:/ECMWFdata/ECME_2020/Apr/ECME_NNL_2020040100/ECME_NNL_202004010000_NL001_LC",sep = "", header = F)

setwd("E:/ECMWFdata/ECME_2020/Nov/ECME_NNL_2020110100")
all_file = list.files()
n = length(all_file)
member = 52
lead_time = 40
file = paste("E:/ECMWFdata/ECME_2020/Nov/ECME_NNL_2020110100", all_file, sep = "/")

df = as.data.frame(matrix(nrow = n*40, ncol = 54))
df2 = as.data.frame(matrix(nrow =242, ncol = 3))

colnames(df)[1] = 'lead-time'
df$`lead-time` = rep(1:40, each = 242)

colnames(df)[2] = 'grid'
df$grid = rep(c(1:242), times = 40)

for (j in 1:n){

  data = read.delim(file[j],sep = "", header = F)
  data_rain = data[c(480:531),c(1:40)]

  data_rain_1 = data_rain[,1:(ncol(data_rain)-1)]
  data_rain_1 = cbind(0,data_rain_1)
  
  data_rain = as.data.frame(lapply(data_rain,as.numeric))
  data_rain_1 = as.data.frame(lapply(data_rain_1,as.numeric))
  
  data_rain = data_rain-data_rain_1
  for (i in 1:40){
  df[(i-1)*242+j,c(3:54)] = data_rain[,i] }

}
for (k in 1:n){
  data = read.delim(file[k],sep = "", header = F)
  
  df2[k,] = data[1,1:3]
  
}




saveRDS(object = mod1, file = "toy model for RR.rds")
mod1 = readRDS(file = "toy model for RR.rds")

r = read.table("D:/ITC/KNMI/radar_data/radar_coord.dat", sep = ",")
r = r[-1,]
r$rID = (1:dim(r)[1])
#rm = c(396*700+369, 396*700+370, 444*700+352, 444*700+353, 434*700+309, 434*700+310, 435*700+309, 435*700+310, 
#       503*700+281, 503*700+282, 503*700+283, 504*700+282, 504*700+283, 448*700+289, 448*700+290, 448*700+291, 
 #      448*700+292, 449*700+288, 449*700+289, 449*700+290, 449*700+291, 449*700+293, 449*700+294, 450*700+288, 
 #      450*700+289, 450*700+294, 451*700+289, 451*700+290, 451*700+291, 451*700+292, 452*700+289, 452*700+290, 
 #      452*700+291, 452*700+292, 325*700+446, 325*700+447, 326*700+446, 326*700+447, 326*700+448, 327*700+446, 
 #      327*700+447, 327*700+448, 329*700+447)

#r = r[-rm,]



r$rLat = r$V2
r$rLon = r$V1
r = r[,3:5]

r = as.data.frame(lapply(r,as.numeric))
df2 = as.data.frame(lapply(df2,as.numeric))
df2$df2id = (1:dim(df2)[1])
df2$df2lat = df2$V2
df2$df2lon = df2$V3
df2 = df2[,4:6]

#d <- pointDistance(df2[,3:2], r[,3:2], lonlat=TRUE, allpairs=T)
d <- pointDistance(r[,3:2], df2[,3:2], lonlat=TRUE, allpairs=T) 
points <- apply(d, 1, which.min)

 library(geosphere)
 distm(c(3.4221, 53.6275), c(3.7643, 53.6275))
 data1 = read.delim("E:/ECMWFdata/ECME_2021/Jan/ECME_NNL_2021010100/ECME_NNL_202101010000_NL001_LC",sep = "", header = F)
 data2 = read.delim("E:/ECMWFdata/ECME_2021/Jan/ECME_NNL_2021010100/ECME_NNL_202101010000_NL002_LC",sep = "", header = F)
 

#df2$ID = r$ID[i]
#df2$distance = d[cbind(1:nrow(d), points)]
#df2$ori_ID = 1:242
r$df2id = points


model.data <- merge(r, df2, by.x="df2id",by.y ="df2id")
model.data = model.data[order(model.data[,2]),]

#setwd("E:/KNMI_DATA")
rvi_table = read.delim("E:/KNMI_DATA/matchdata.txt", sep = "\t")

library(raster)
library(ncdf4)

ra = raster("D:/ITC/KNMI/radar_data/09-20211013T195924Z-001/09/RAD_NL25_RAC_MFBS_01H_202009010200.nc")
ra = as.data.frame(ra)
ra$raid = c(1:535500)
ra_nl = na.omit(ra)
nl.data <- merge(model.data, ra_nl, by.x="rID",by.y ="raid")

setwd("E:/KNMI_DATA")
write.table(model.data, "matchdata.txt", sep = "\t", col.names = T, row.names = F)
#####Basic function for working on R#########################################################
#1. 批量更改文件夹里的文件名
setwd("E:/ECMWFdata/pre-process/Y2021")
files<-list.files()

for (f in files){
  newname<-sub('Apr','04',f) #将原文件中的字符xt，替换为字符ab
  file.rename(f,newname)
}
#2. 忽略NA
na.rm = T
#3. ggplot作图

ggplot(r %>% mutate(celln = 1:nrow(r))) + geom_point(aes(x = V3, y = V2, col = celln))
head(as.data.frame(ra, xy=TRUE))

p <- ggplot(ALL_5_3,aes(x=x, y=ALL_5_3_500_bs))+geom_point(size = 1, shape = 21)+xlab("p0")+ylab("Observation")+
    xlim(0,15) + ylim(0,0.5) + ggtitle("Observation versus mean for 60 hours forecast time (cube root)") 
p <- ggplot(final,aes(x=EC4meanfinal1, y=R4final1))+geom_point(size = 1, shape = 21)+xlab("p0")+ylab("Observation")+
  xlim(0, 20) + ylim(0,20) + ggtitle("Observation versus mean for 60 hours forecast time") 


set.seed(1234)
threshold <-  c(0.05,0.5,1,3,5,10,15)
type <- rep(c('3var_bs','5var_bs','default_bs'),each = 7)

value <- c(0.1025,0.0782,0.0661,0.0356,0.0199,0.00389,0.000279,0.1028,0.0782
           ,0.0662,0.0356,0.0198,0.00392,0.000278,0.1058,0.0798,0.0677,0.03600,0.02002,0.003858,0.0002769)
df <- data.frame(lead_time = lead_time, type = type, value = value)
ggplot(data = df, mapping = aes(x = lead_time, y = value, colour = type)) + geom_line()+xlab("lead time (h)")+ylab("crps")+
  ggtitle("Continuous ranked probability score of emos (ZAGA) and qrf with different lead time") 


type <- rep(c('3var_ss','5var_ss','default_ss'),each = 7)
value <- c(0.499,0.474,0.445,0.366,0.267,0.0983,0.003901,0.498,0.474,0.444,0.365,0.2678,0.0905,0.008611,0.483,
           0.463,0.432,0.358,0.261,0.1056,0.01185)

df <- data.frame(threshold = threshold, type = type, value = value)
ggplot(data = df, mapping = aes(x = lead_time, y = value, colour = type)) + geom_bar(stat = 'identity')



x <- rep(c(24,48,72,96), each = 3)
y <- rep(c('raw','zaga','qrf'),times = 4)
set.seed(1234)
bs <- c(0.1025,0.1028,0.1058,0.07819895,0.07822975,0.07981754,0.06609,0.0662467,0.0677281,0.035556,0.035612,
       0.036003,0.019851,0.01983735,0.02002396,0.003889,0.00392298,0.003857852,0.000279,0.00027778,0.0002769)
#z <- c()
df <- data.frame(x = x, y=y, value1=value1)
#不作任何条形宽度和条形距离的调整
ggplot(data = df, mapping = aes(x = factor(x), y = value1, fill = y)) + geom_bar(stat = 'identity', position = 'dodge')+xlab("lead time (h)")+ylab("crps")+
  ggtitle("Continuous ranked probability score of emos (ZAGA) and qrf with different lead time") 


ss <- c(0.4990924,0.4976741,0.4826699,0.4742184,0.4740113,0.4633356,0.4453593,0.4440599,0.4316278,0.3656134,
        0.3646171,0.3578695,0.2672695,0.2677825,0.2608949,0.09827222,0.09048656,0.1055856,0.003901127,0.008611494,0.01184584)
df <- data.frame(x = x, y=y,ss=ss)
#不作任何条形宽度和条形距离的调整
ggplot(data = df, mapping = aes(x = factor(x), y = ss, fill = y)) + geom_bar(stat = 'identity', position = 'dodge')+xlab("threshold")+
  ggtitle("Comparasion of ss for qrf (48h)") 



#4.写入/读取mod格式文件
saveRDS(object = mod1, file = "24h_forecast_3rd_year_2step.rds")
setwd("E:/output")
mod1 = readRDS(file = "24h_forecast_3rd_year_3step.rds")
####check with radar data-------------------------------------------------------------------------------






setwd("E:/RadarData/2019/07")

all_file = list.files()
n = length(all_file)
file = paste("E:/RadarData/2019/07", all_file, sep = "/")

df = as.data.frame(matrix(nrow = n, ncol = 1))
#df = as.data.frame(matrix(nrow = n, ncol =38063))

for (i in 1:n){
 ra = raster(file[i])
 ra = as.data.frame(ra)
 ra = na.omit(ra)
 max = max(ra)
# rain = ra$image1_image_data
# df[i,]=rain
  df[i,]= max
 }
max(df)


  




####ECMWF_data_organize---------------------------------------------------------------------------------------------


#initial time, grid, lead-time, member1 to member 50

  df2[c(((i-1)*242+1):(i*242)),] = b

colnames(b)[53] = 'lead-time'
b$`lead-time` = 1:40
colnames(b)[54] = 'grid'
b$grid = rep(c(1:242), each = 40)


######working with nc_data--------------------------------------------------------------------------------------------------
#working with date time





library(ncdf4)
library(raster)
rb = nc_open("D:/ITC/KNMI/radar_data/09-20211013T195924Z-001/09/RAD_NL25_RAC_MFBS_01H_202009010200.nc")
ra = raster("E:/RadarData/2021/02/RAD_NL25_RAC_MFBS_01H_202102030700.nc")+
  raster("E:/RadarData/2021/02/RAD_NL25_RAC_MFBS_01H_202102030800.nc")+
  raster("E:/RadarData/2021/02/RAD_NL25_RAC_MFBS_01H_202102030900.nc")+
  raster("E:/RadarData/2021/02/RAD_NL25_RAC_MFBS_01H_202102031000.nc")+
  raster("E:/RadarData/2021/02/RAD_NL25_RAC_MFBS_01H_202102031100.nc")+
  raster("E:/RadarData/2021/02/RAD_NL25_RAC_MFBS_01H_202102031200.nc")

plot(ra, main = "Radar Precipitation for 6 hours (mm)", xlim = c(200,600), ylim = c(-4300,-3800))





ra = as.data.frame(ra)
ra$ID = (1:dim(ra)[1])
ra = na.omit(ra)

r <- raster(nrow=4, ncol=5)
values(r) <- 1:ncell(r)
s <- raster(nrow=255, ncol=233)
s <- resample(ra, s, method='ngb')

matchdata <- merge(model.data, ra, by.x="rID",by.y ="ID")

nc_data = nc_open("D:/ITC/thesis/from_10_29/SACA&D_data/rr_0.25deg_reg_2011-2017_v2.0_saobs.nc/rr_0.25deg_reg_2011-2017_v2.0_saobs.nc")
ra = nc_open("D:/ITC/KNMI/radar_data/09-20211013T195924Z-001/09/RAD_NL25_RAC_MFBS_01H_202009010200.nc")
{ 
  sink ('test.txt')
  print(ra)
  sink()
  }

rb = raster("D:/ITC/KNMI/radar_data/09-20211013T195924Z-001/09/RAD_NL25_RAC_MFBS_01H_202009010200.nc")

ra = nc_open("E:/RadarData/RADNL_CLIM____MFBSNL25_01H_20190101T000000_20200101T000000_netCDF4_0002/2019/11/RAD_NL25_RAC_MFBS_01H_201911010300.nc")
rc = nc_open("D:/ITC/KNMI/radar_data/09-20211013T195924Z-001/09/RAD_NL25_RAC_MFBS_01H_202009070200.nc")


lon <- ncvar_get(ra, "x")
lat <- ncvar_get(ra, "y", verbose = F)

try <- ncvar_get(ra,var)

t <- ncvar_get(rb, "time")
tunits <- ncatt_get(rb,"time","units")
attributes(nc$var)$names
ncatt_get(rb, attributes(rb$var)$names[3])
data1 <- ncvar_get(rc, attributes(rc$var)$names[3])
ncatt_get(rc, attributes(ra$var)$names[3])


rr.array <- ncvar_get(ra, "projection")
fillvalue <- ncatt_get(nc_data, "", "_FillValue")
fillvalue
nc_close(nc_data) 
rr.array[rr.array == fillvalue$value] <- NA
rr.slice <- rr.array[, , 100] 
dim(rr.slice)
r <- raster(t(rr.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
crs(r)
r <- flip(r, direction='y')
plot(r)

lat = lat[67:77,]
lon = lon[102:139,]






####reading radar data and reprojection----------------------------------------------------------------------------------------------------------------------------


library(rgdal)

library(ncdf4)


#nc <- nc_open("D:/ITC/KNMI/radar_data/09-20211013T195924Z-001/09/RAD_NL25_RAC_MFBS_01H_202009010200.nc")
ra = raster("D:/ITC/KNMI/radar_data/09-20211013T195924Z-001/09/RAD_NL25_RAC_MFBS_01H_202009010200.nc")
ra2 = raster("D:/ITC/KNMI/radar_data/09-20211013T195924Z-001/09/RAD_NL25_RAC_MFBS_01H_202009151200.nc")
ra3 = raster("D:/ITC/KNMI/radar_data/09-20211013T195924Z-001/09/RAD_NL25_RAC_MFBS_01H_202009270800.nc")
r = read.table("D:/ITC/KNMI/radar_data/radar_coord.dat", sep = ",")
r = r[(2:535501),]

setwd("E:/ECMWFdata/ECME_NNL_2020110100")
all_file = list.files()
n = length(all_file)

df2 = as.data.frame(matrix(nrow =242, ncol = 3))


for (k in 1:n){
  data = read.delim(file[k],sep = "", header = F)

  df2[k,] = data[1,1:3]
  
}




proj4string(ra) <- "+proj=stere +lat_0=90 +lon_0=0 +k=0.933012709177451 +units=km +datum=WGS84"
ra = projectRaster(from=ra, crs=CRS("+init=epsg:4326"), method = "ngb")

proj4string(ra2) <- "+proj=stere +lat_0=90 +lon_0=0 +k=0.933012709177451 +units=km +datum=WGS84"
ra2 = projectRaster(from=ra2, crs=CRS("+init=epsg:4326"), method = "ngb")

proj4string(ra) <- "+proj=stere +lat_0=90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6356752 +units=m +no_defs"
ra = projectRaster(from=ra, crs=CRS("+init=epsg:4326"), method = "ngb")

crs = "+proj=stere +lat_0=90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6356752 +units=m +no_defs "
a = projectRaster(ra,crs=crs)

b = as.data.frame(ra, xy = TRUE)
b = na.omit(b)
b2 = as.data.frame(ra2,xy=TRUE)

plot(ra)
crs(ra)
#crs1="+proj=utm +zone=32N +datum=WGS84 +units=m +no_defs"
#crs2 = "+proj=utm +zone=19 +datum=NAD83 +units=m +no_defs"
#crs3 = "+ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
crs3 = "+ellps=WGS84 +datum=WGS84 +no_defs+towgs84=0,0,0"
#crs4 = "+proj=longlat +datum=WGS84 +no_defs"

a=projectRaster(ra, crs =crs3,rm.omit=F)

library(rgdal)

r.pts <- rasterToPoints(ra, spatial=TRUE)
proj4string(r.pts)

projection(x) <- crs("+ellps=WGS84 +datum=WGS84 +no_defs+towgs84=0,0,0")

plot(a)
crs(a)
crs(ra)

getwd()
setwd(E:/ECMWFdata/ECME_NNL_2018110100/)

#initial time, grid, lead-time, member1 to member 50
40*50*242*

  

ggplot(r %>% mutate(celln = 1:nrow(r))) + geom_point(aes(x = V3, y = V2, col = celln))
head(as.data.frame(ra, xy=TRUE))

####resample the radar data #########################################
ra = raster("E:/RadarData/2019/01/RAD_NL25_RAC_MFBS_01H_201901010100.nc")
#remove clutter points
ra[397,369] = ra[397,370]= ra[445,352] = ra[445,353] = ra[435,309] = ra[435,310] =
  ra[436,309] = ra[436,310] = ra[504,281] = ra[504,282] = ra[504,283] = ra[505,282] =
  ra[505,283] = ra[449,289] = ra[449,290] = ra[449,291] = ra[449,292] = ra[450,288] = 
  ra[450,289] = ra[450,290] = ra[450,291] = ra[450,293] = ra[450,294] = ra[451,288] = 
  ra[451,289] = ra[451,294] = ra[452,289] = ra[452,290] = ra[452,291] = ra[452,292] = 
  ra[453,289] = ra[453,290] = ra[453,291] = ra[453,292] = ra[326,446] = ra[326,447] = 
  ra[327,446] = ra[327,447] = ra[327,448] = ra[328,446] = ra[328,447] = ra[328,448] = 
  ra[330,447] = NA


#ra[is.na(ra[])] <- 0  #replace all NA values in a raster by 0
lat <- init(ra, 'y')
lat = as.data.frame(lat)
colnames(lat)[1] <- 'lat'
lon <- init(ra, 'x')
lon = as.data.frame(lon)
colnames(lon)[1] <- 'lon'
location = cbind(lat,lon)
location$ID = (1:dim(location)[1])
#new_ra1 = aggregate(ra,fact=3)
#max RA 9*9
new_ra1 = aggregate(ra, fact=9, expand=FALSE, fun=mean, na.rm=TRUE)
#new_ra1 = aggregate(new_ra, fact=3, expand=FALSE, fun=max, na.rm=TRUE)

lat = init(ra,'y')
lat = as.data.frame(lat)



data_matrix <- rasterToPoints(new_ra1)
data_matrix= as.data.frame(data_matrix)
colnames(data_matrix)[1] <- 'lon'
colnames(data_matrix)[2] <- 'lat'
new_ra = merge(data_matrix, location, by = c("lat","lon")) 
rvi_table = read.delim("E:/KNMI_DATA/matchdata.txt", sep = "\t")
radar_output = merge(new_ra,rvi_table,by.x = "ID", by.y = "rID")

#Make the list of ID for matching
upscal_match_9km = cbind.data.frame(radar_output$ID,radar_output$df2id)
colnames(upscal_match_9km) = c("rID","dfID")
#write.table(upscal_match_9km, "upscal_match_9km.txt", sep = "\t", col.names = T, row.names = F)

ggplot(radar_output) +
  geom_point(aes(x = rLon, y= rLat, color = image1_image_data)) 





setwd("E:/RadarData/2018")
first_category_name = list.files("12")
dir = paste("E:/RadarData/2018/12/",first_category_name,sep="")
n = length(dir) 

radar_combine = as.data.frame(matrix(0,nrow = 599, ncol = n))
#date = as.data.frame(matrix(0,nrow = 1, ncol = n))
#names(radar_combine) = c(1:n)

for(i in 1:n){
  ra<-raster(dir[i]) 
#  rb<-nc_open(dir[i])
  ra[397,369] = ra[397,370]= ra[445,352] = ra[445,353] = ra[435,309] = ra[435,310] =
    ra[436,309] = ra[436,310] = ra[504,281] = ra[504,282] = ra[504,283] = ra[505,282] =
    ra[505,283] = ra[449,289] = ra[449,290] = ra[449,291] = ra[449,292] = ra[450,288] = 
    ra[450,289] = ra[450,290] = ra[450,291] = ra[450,293] = ra[450,294] = ra[451,288] = 
    ra[451,289] = ra[451,294] = ra[452,289] = ra[452,290] = ra[452,291] = ra[452,292] = 
    ra[453,289] = ra[453,290] = ra[453,291] = ra[453,292] = ra[326,446] = ra[326,447] = 
    ra[327,446] = ra[327,447] = ra[327,448] = ra[328,446] = ra[328,447] = ra[328,448] = 
    ra[330,447] = NA
  
 #get 3*3 
  new_ra1 = aggregate(ra, fact=9, expand=FALSE, fun=mean, na.rm=TRUE)
  #get the resolution same as ECMWF file
  #aggregate new_ra to the resolution of ECMWF resolution , fun=max
  
 # new_ra1 = aggregate(new_ra, fact=3, expand=FALSE, fun=max, na.rm=TRUE)
  
  data_matrix <- rasterToPoints(new_ra1)
  data_matrix= as.data.frame(data_matrix)
  output = data_matrix$image1_image_data
  radar_combine[,i]=output
#  dates <- as.POSIXct(ncvar_get(rb, "time"), origin = "2000-01-01 01:00", tz = "UTC")
#  date[,i]=dates
}



new_df = radar_combine

f <- matrix(0, dim(new_df)[1], dim(new_df)[2] %/% 6 )
for (j in seq(1, (dim(new_df)[2] - dim(new_df)[2] %% 6), by = 6)) {
  f[, (j - 1) / 6 + 1] <- rowSums(new_df[, j:(j + 5)],na.rm=TRUE)
}

f = as.data.frame(f)
colnames(f) = 1: dim(f)[2]
setwd("E:/output/radar_output")

write.table(f, "2020_11_9km_avg.txt", sep = "\t", col.names = T, row.names = F)




#############################################################################################
#pre-processing data
#####(for WDIR)######################################################################################

setwd('E:/ECMWFdata/ECME_2020')
# data = read.delim("E:/ECMWFdata/ECME_2020/Jun/ECME_NNL_2020062112/ECME_NNL_202006211200_NL224_LC",sep = "", header = F)
first_category_name = list.files("Oct")
dir = paste("./Oct/",first_category_name,sep="")
n = length(dir) 
n_sub<-rep(0,n)
n_sub<-as.data.frame(n_sub)
n_sub<-t(n_sub)
head(n_sub)
#create the df for HRES
merge_HRES = as.data.frame(matrix(1:40,nrow = 1, ncol = 40))
names(merge_HRES) = c(1:40)
merge_HRES$second_category <- 'second_category'
merge_HRES$first_category <- 'first_category'
#create the df for avg_sin
merge_avg_sin = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_avg_sin) = c(1:60)
merge_avg_sin$second_category <- 'second_category'
merge_avg_sin$first_category <- 'first_category'
#create the df for avg_cos
merge_avg_cos = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_avg_cos) = c(1:60)
merge_avg_cos$second_category <- 'second_category'
merge_avg_cos$first_category <- 'first_category'
#create the df for median_sin
merge_median_sin = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_median_sin) = c(1:60)
merge_median_sin$second_category <- 'second_category'
merge_median_sin$first_category <- 'first_category'
#create the df for median_cos
merge_median_cos = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_median_cos) = c(1:60)
merge_median_cos$second_category <- 'second_category'
merge_median_cos$first_category <- 'first_category'
#create the df for q10_sin
merge_q10_sin = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q10_sin) = c(1:60)
merge_q10_sin$second_category <- 'second_category'
merge_q10_sin$first_category <- 'first_category'
#create the df for q10_cos
merge_q10_cos = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q10_cos) = c(1:60)
merge_q10_cos$second_category <- 'second_category'
merge_q10_cos$first_category <- 'first_category'
#create the df for q90_sin
merge_q90_sin = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q90_sin) = c(1:60)
merge_q90_sin$second_category <- 'second_category'
merge_q90_sin$first_category <- 'first_category'
#create the df for q90_cos
merge_q90_cos = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q90_cos) = c(1:60)
merge_q90_cos$second_category <- 'second_category'
merge_q90_cos$first_category <- 'first_category'



for(i in 1:n){         #对于每个一级目录(文件夹)
  b=list.files(dir[i]) #b是列出每个一级目录(文件夹)中每个xlsx文件的名称
  n_sub[i]=length(b)   #得到一级目录(文件夹)下xlsx的文件个数:n_sub
  
  for(j in 1:n_sub[i]){     #对于每个一级目录(文件夹)下的每个xlsx文件
    data<-read.delim(file=paste(dir[i],'/',b[j],sep=''),sep = "", header = F) #读取xlsx文件
    #WDIR
    data_WDIR = data[108:160,]
    data_WDIR[,48:83]=trunc(data_WDIR[,48:83]/1000000)
    data_WDIR[2,64]=data_WDIR[2,64]/10^234
    WDIR_6hstep = data_WDIR[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
    WDIR_HRES = WDIR_6hstep[2,1:40]
    WDIR_ENS = WDIR_6hstep[3:53,]
    WDIR_ENS = as.data.frame(lapply(WDIR_ENS,as.numeric))
    avg_var_sin_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    avg_var_cos_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    median_var_sin_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    median_var_cos_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q10_var_sin_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q10_var_cos_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q90_var_sin_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q90_var_cos_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    for(k in 1:ncol(WDIR_ENS)){ 
      avg_var <- mean(WDIR_ENS[, k])
      avg_var_sin <- sin(avg_var*pi/180)
      avg_var_cos <- cos(avg_var*pi/180)
      median_var <- median(WDIR_ENS[,k])
      median_var_sin <- sin(median_var*pi/180)
      median_var_cos <- cos(median_var*pi/180)
      q10_var = quantile(WDIR_ENS[,k],probs=0.1)
      q10_var_sin <- sin(q10_var*pi/180)
      q10_var_cos <- cos(q10_var*pi/180)
      q90_var = quantile(WDIR_ENS[,k],probs=0.9)
      q90_var_sin <- sin(q90_var*pi/180)
      q90_var_cos <- cos(q90_var*pi/180)
      avg_var_sin_df[,k]=avg_var_sin
      avg_var_cos_df[,k]=avg_var_cos
      median_var_sin_df[,k]=median_var_sin
      median_var_cos_df[,k]=median_var_cos
      q10_var_sin_df[,k]=q10_var_sin
      q10_var_cos_df[,k]=q10_var_cos
      q90_var_sin_df[,k]=q90_var_sin
      q90_var_cos_df[,k]=q90_var_cos
    }    
    #new row for HRES
    new_HRES = WDIR_HRES
    names(new_HRES)<-c(1:40)
     new_HRES$second_category<-substr(b[j],14,27)      
    new_HRES$first_category<-first_category_name[i]   
    merge_HRES<-rbind(merge_HRES,new_HRES)
    #new row for avg_sin
    new_avg_sin = avg_var_sin_df
    names(new_avg_sin)<-c(1:60)
    new_avg_sin$second_category<-substr(b[j],14,27)       
    new_avg_sin$first_category<-first_category_name[i] 
    merge_avg_sin<-rbind(merge_avg_sin,new_avg_sin)
    #new row for avg_cos
    new_avg_cos = avg_var_cos_df
    names(new_avg_cos)<-c(1:60)
    new_avg_cos$second_category<-substr(b[j],14,27)       
    new_avg_cos$first_category<-first_category_name[i] 
    merge_avg_cos<-rbind(merge_avg_cos,new_avg_cos)    
    #new row for median_sin
    new_median_sin = median_var_sin_df
    names(new_median_sin)<-c(1:60)
    new_median_sin$second_category<-substr(b[j],14,27)       
    new_median_sin$first_category<-first_category_name[i] 
    merge_median_sin<-rbind(merge_median_sin,new_median_sin)
    #new row for median_cos
    new_median_cos = median_var_cos_df
    names(new_median_cos)<-c(1:60)
    new_median_cos$second_category<-substr(b[j],14,27)       
    new_median_cos$first_category<-first_category_name[i] 
    merge_median_cos<-rbind(merge_median_cos,new_median_cos)
    #new row for q10_sin
    new_q10_sin = q10_var_sin_df
    names(new_q10_sin)<-c(1:60)
    new_q10_sin$second_category<-substr(b[j],14,27)       
    new_q10_sin$first_category<-first_category_name[i] 
    merge_q10_sin<-rbind(merge_q10_sin,new_q10_sin)
    #new row for q10_cos
    new_q10_cos = q10_var_cos_df
    names(new_q10_cos)<-c(1:60)
    new_q10_cos$second_category<-substr(b[j],14,27)       
    new_q10_cos$first_category<-first_category_name[i] 
    merge_q10_cos<-rbind(merge_q10_cos,new_q10_cos)
    #new row for q90_sin
    new_q90_sin = q90_var_sin_df
    names(new_q90_sin)<-c(1:60)
    new_q90_sin$second_category<-substr(b[j],14,27)       
    new_q90_sin$first_category<-first_category_name[i] 
    merge_q90_sin<-rbind(merge_q90_sin,new_q90_sin)
    #new row for q90_cos
    new_q90_cos = q90_var_cos_df
    names(new_q90_cos)<-c(1:60)
    new_q90_cos$second_category<-substr(b[j],14,27)       
    new_q90_cos$first_category<-first_category_name[i] 
    merge_q90_cos<-rbind(merge_q90_cos,new_q90_cos)
  }
}       
setwd("E:/ECMWFdata/pre-process")
write.table(merge_HRES, "WDIR_6h_HRES.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_avg_sin, "WDIR_6h_mean_sin.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_avg_cos, "WDIR_6h_mean_cos.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_median_sin, "WDIR_6h_median_sin.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_median_cos, "WDIR_6h_median_cos.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q10_sin, "WDIR_6h_q10_sin.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q10_cos, "WDIR_6h_q10_cos.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q90_sin, "WDIR_6h_q90_sin.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q90_cos, "WDIR_6h_q90_cos.txt", sep = "\t", col.names = T, row.names = F)



####(for_RR) (mm)###############################################################################################    
library(e1071)
setwd('E:/ECMWFdata/ECME_2020')
# data = read.delim("E:/ECMWFdata/ECME_2019/Dec/ECME_NNL_2019121512/ECME_NNL_201912151200_NL084_LC",sep = "", header = F)
first_category_name = list.files("Feb")
dir = paste("./Feb/",first_category_name,sep="")
n = length(dir) 
n_sub<-rep(0,n)
n_sub<-as.data.frame(n_sub)
n_sub<-t(n_sub)
head(n_sub)

#create the df for HRES
merge_HRES = as.data.frame(matrix(1:40,nrow = 1, ncol = 40))
names(merge_HRES) = c(1:40)
merge_HRES$second_category <- 'second_category'
merge_HRES$first_category <- 'first_category'
#create the df for avg
merge_avg = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_avg) = c(1:60)
merge_avg$second_category <- 'second_category'
merge_avg$first_category <- 'first_category'
#create the df for sd
merge_sd = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_sd) = c(1:60)
merge_sd$second_category <- 'second_category'
merge_sd$first_category <- 'first_category'
#create the df for median
merge_median = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_median) = c(1:60)
merge_median$second_category <- 'second_category'
merge_median$first_category <- 'first_category'
#create the df for q10
merge_q10 = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q10) = c(1:60)
merge_q10$second_category <- 'second_category'
merge_q10$first_category <- 'first_category'
#create the df for q90
merge_q90 = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q90) = c(1:60)
merge_q90$second_category <- 'second_category'
merge_q90$first_category <- 'first_category'
#create the df for p0
merge_p0 = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_p0) = c(1:60)
merge_p0$second_category <- 'second_category'
merge_p0$first_category <- 'first_category'
#create the df for p1
merge_p1 = as.data.frame(matrix(1:60, nrow = 1, ncol = 60))
names(merge_p1) = c(1:60)
merge_p1$second_category <- 'second_category'
merge_p1$first_category <- 'first_category'
#create the df for p3
merge_p3 = as.data.frame(matrix(1:60, nrow = 1, ncol = 60))
names(merge_p3) = c(1:60)
merge_p3$second_category <- 'second_category'
merge_p3$first_category <- 'first_category'
#create the df for p5
merge_p5 = as.data.frame(matrix(1:60, nrow = 1, ncol = 60))
names(merge_p5) = c(1:60)
merge_p5$second_category <- 'second_category'
merge_p5$first_category <- 'first_category'
#create the df for p10
merge_p10 = as.data.frame(matrix(1:60, nrow = 1, ncol = 60))
names(merge_p10) = c(1:60)
merge_p10$second_category <- 'second_category'
merge_p10$first_category <- 'first_category'
#create the df for p20
merge_p20 = as.data.frame(matrix(1:60, nrow = 1, ncol = 60))
names(merge_p20) = c(1:60)
merge_p20$second_category <- 'second_category'
merge_p20$first_category <- 'first_category'
#create the df for skewness
merge_skewness = as.data.frame(matrix(1:60, nrow = 1, ncol = 60))
names(merge_skewness) = c(1:60)
merge_skewness$second_category <- 'second_category'
merge_skewness$first_category <- 'first_category'
#create the df for kurtosis
merge_kurtosis = as.data.frame(matrix(1:60, nrow = 1, ncol = 60))
names(merge_kurtosis) = c(1:60)
merge_kurtosis$second_category <- 'second_category'
merge_kurtosis$first_category <- 'first_category'


for(i in 1:n){         #对于每个一级目录(文件夹)
  b=list.files(dir[i]) #b是列出每个一级目录(文件夹)中每个xlsx文件的名称
  n_sub[i]=length(b)   #得到一级目录(文件夹)下xlsx的文件个数:n_sub
  
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
    RR_HRES = RR_6hstep[2,1:40]
    RR_ENS = RR_6hstep[3:53,]
    RR_ENS = as.data.frame(lapply(RR_ENS,as.numeric))
    avg_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    sd_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    median_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q10_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q90_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    p0_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    p1_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    p3_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    p5_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    p10_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    p20_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    skewness_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    kurtosis_var_df <- as.data.frame(matrix(1:60, nrow = 1, ncol = 60))
    for(k in 1:ncol(RR_ENS)){ 
      #  p0 <- colSums(RR_ENS>0)/length(RR_ENS[,k])   
      avg_var_df[,k] <- mean(RR_ENS[, k])
      sd_var_df[,k] <- sd(RR_ENS[,k])
      median_var_df[,k] <- median(RR_ENS[,k])
      q10_var_df[,k] = quantile(RR_ENS[,k],probs=0.1)
#      0.25
#      0.75

      q90_var_df[,k] = quantile(RR_ENS[,k],probs=0.9)
      p0_var_df[,k] = sum(RR_ENS[,k] > 0)/length(RR_ENS[,k])
      p1_var_df[,k] = sum(RR_ENS[,k] > 1)/length(RR_ENS[,k])
      p3_var_df[,k] = sum(RR_ENS[,k] > 3)/length(RR_ENS[,k])
      p5_var_df[,k] = sum(RR_ENS[,k] > 5)/length(RR_ENS[,k])
      p10_var_df[,k] = sum(RR_ENS[,k] > 10)/length(RR_ENS[,k])
      p20_var_df[,k] = sum(RR_ENS[,k] > 20)/length(RR_ENS[,k])
      skewness_var_df[,k] = skewness(RR_ENS[,k])
      kurtosis_var_df[,k] = kurtosis(RR_ENS[,k])
    }    
    
    
    #new row for HRES
#    new_HRES = RR_HRES
    names(RR_HRES)<-c(1:40)
    RR_HRES$second_category<-substr(b[j],14,27)      
    RR_HRES$first_category<-first_category_name[i]   
    merge_HRES<-rbind(merge_HRES,RR_HRES)
    #new row for avg
#    new_avg = avg_var_df
    names(avg_var_df)<-c(1:60)
    avg_var_df$second_category<-substr(b[j],14,27)       
    avg_var_df$first_category<-first_category_name[i] 
    merge_avg<-rbind(merge_avg,avg_var_df)    
    #new row for sd
#    new_sd = sd_var_df
    names(sd_var_df)<-c(1:60)
    sd_var_df$second_category<-substr(b[j],14,27)       
    sd_var_df$first_category<-first_category_name[i] 
    merge_sd<-rbind(merge_sd,sd_var_df) 
    #new row for median
#    new_median = median_var_df
    names(median_var_df)<-c(1:60)
    median_var_df$second_category<-substr(b[j],14,27)       
    median_var_df$first_category<-first_category_name[i] 
    merge_median<-rbind(merge_median,median_var_df)
    #new row for q10
#    new_q10 = q10_var_df
    names(q10_var_df)<-c(1:60)
    q10_var_df$second_category<-substr(b[j],14,27)       
    q10_var_df$first_category<-first_category_name[i] 
    merge_q10<-rbind(merge_q10,q10_var_df)
    #new row for q90
#    new_q90 = q90_var_df
    names(q90_var_df)<-c(1:60)
    q90_var_df$second_category<-substr(b[j],14,27)       
    q90_var_df$first_category<-first_category_name[i] 
    merge_q90<-rbind(merge_q90,q90_var_df)
    #new row for p0
#    new_p0 = p0_var_df
    names(p0_var_df)<-c(1:60)
    p0_var_df$second_category<-substr(b[j],14,27)       
    p0_var_df$first_category<-first_category_name[i] 
    merge_p0<-rbind(merge_p0,p0_var_df)    
    #new row for p1
#    new_p1 = p1_var_df
    names(p1_var_df)<-c(1:60)
    p1_var_df$second_category<-substr(b[j],14,27)       
    p1_var_df$first_category<-first_category_name[i] 
    merge_p1<-rbind(merge_p1,p1_var_df)   
    #new row for p3
#    new_p3 = p3_var_df
    names(p3_var_df)<-c(1:60)
    p3_var_df$second_category<-substr(b[j],14,27)       
    p3_var_df$first_category<-first_category_name[i] 
    merge_p3<-rbind(merge_p3,p3_var_df)  
    #new row for p5
#    new_p5 = p5_var_df
    names(p5_var_df)<-c(1:60)
    p5_var_df$second_category<-substr(b[j],14,27)       
    p5_var_df$first_category<-first_category_name[i] 
    merge_p5<-rbind(merge_p5,p5_var_df) 
    #new row for p10
#    new_p10 = p10_var_df
    names(p10_var_df)<-c(1:60)
    p10_var_df$second_category<-substr(b[j],14,27)       
    p10_var_df$first_category<-first_category_name[i] 
    merge_p10<-rbind(merge_p10,p10_var_df)  
    #new row for p20
#    new_p20 = p20_var_df
    names(p20_var_df)<-c(1:60)
    p20_var_df$second_category<-substr(b[j],14,27)       
    p20_var_df$first_category<-first_category_name[i] 
    merge_p20<-rbind(merge_p20,p20_var_df) 
    #new row for skewness
#    new_skewness = skewness_var_df
    names(skewness_var_df)<-c(1:60)
    skewness_var_df$second_category<-substr(b[j],14,27)       
    skewness_var_df$first_category<-first_category_name[i] 
    merge_skewness<-rbind(merge_skewness,skewness_var_df)    
    #new row for kurtosis
#    new_kurtosis = kurtosis_var_df
    names(kurtosis_var_df)<-c(1:60)
    kurtosis_var_df$second_category<-substr(b[j],14,27)       
    kurtosis_var_df$first_category<-first_category_name[i] 
    merge_kurtosis<-rbind(merge_kurtosis,kurtosis_var_df)    
  }
}    

setwd("E:/ECMWFdata/pre-process")
write.table(merge_HRES, "RR_6h_HRES_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_avg, "RR_6h_mean_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_sd, "RR_6h_sd_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_median, "RR_6h_median_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q10, "RR_6h_q10_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q90, "RR_6h_q90_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_p0, "RR_6h_p0_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_p1, "RR_6h_p1_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_p3, "RR_6h_p3_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_p5, "RR_6h_p5_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_p10, "RR_6h_p10_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_p20, "RR_6h_p20_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_skewness, "RR_6h_skewness_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_kurtosis, "RR_6h_kurtosis_2018.txt", sep = "\t", col.names = T, row.names = F)

#####（for_convectiveRR)##########################################################################
library(e1071)
Sys.time()
setwd('E:/ECMWFdata/ECME_2018')
# data = read.delim("E:/ECMWFdata/ECME_2019/Dec/ECME_NNL_2019121512/ECME_NNL_201912151200_NL084_LC",sep = "", header = F)
first_category_name = list.files("Dec")
dir = paste("./Dec/",first_category_name,sep="")
n = length(dir) 
n_sub<-rep(0,n)
n_sub<-as.data.frame(n_sub)
n_sub<-t(n_sub)
head(n_sub)
# create df2 for conRR
#create the df for HRES
merge_HRES_conRR = as.data.frame(matrix(1:40,nrow = 1, ncol = 40))
names(merge_HRES_conRR) = c(1:40)
merge_HRES_conRR$second_category <- 'second_category'
merge_HRES_conRR$first_category <- 'first_category'
#create the df for avg
merge_avg_conRR = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_avg_conRR) = c(1:60)
merge_avg_conRR$second_category <- 'second_category'
merge_avg_conRR$first_category <- 'first_category'
#create the df for sd
merge_sd_conRR = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_sd_conRR) = c(1:60)
merge_sd_conRR$second_category <- 'second_category'
merge_sd_conRR$first_category <- 'first_category'
#create the df for median
merge_median_conRR = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_median_conRR) = c(1:60)
merge_median_conRR$second_category <- 'second_category'
merge_median_conRR$first_category <- 'first_category'
#create the df for q10
merge_q10_conRR = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q10_conRR) = c(1:60)
merge_q10_conRR$second_category <- 'second_category'
merge_q10_conRR$first_category <- 'first_category'
#create the df for q90
merge_q90_conRR = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q90_conRR) = c(1:60)
merge_q90_conRR$second_category <- 'second_category'
merge_q90_conRR$first_category <- 'first_category'
#create the df for q25
merge_q25_conRR = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q25_conRR) = c(1:60)
merge_q25_conRR$second_category <- 'second_category'
merge_q25_conRR$first_category <- 'first_category'
#create the df for q75
merge_q75_conRR = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q75_conRR) = c(1:60)
merge_q75_conRR$second_category <- 'second_category'
merge_q75_conRR$first_category <- 'first_category'
#create the df for p0.05
merge_p0.05_conRR = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_p0.05_conRR) = c(1:60)
merge_p0.05_conRR$second_category <- 'second_category'
merge_p0.05_conRR$first_category <- 'first_category'
#create the df for p0.3
merge_p0.3_conRR = as.data.frame(matrix(1:60, nrow = 1, ncol = 60))
names(merge_p0.3_conRR) = c(1:60)
merge_p0.3_conRR$second_category <- 'second_category'
merge_p0.3_conRR$first_category <- 'first_category'
#create the df for p10
merge_p10_conRR = as.data.frame(matrix(1:60, nrow = 1, ncol = 60))
names(merge_p10_conRR) = c(1:60)
merge_p10_conRR$second_category <- 'second_category'
merge_p10_conRR$first_category <- 'first_category'

#create the df for min
merge_min_conRR = as.data.frame(matrix(1:60, nrow = 1, ncol = 60))
names(merge_min_conRR) = c(1:60)
merge_min_conRR$second_category <- 'second_category'
merge_min_conRR$first_category <- 'first_category'

#create the df for max
merge_max_conRR = as.data.frame(matrix(1:60, nrow = 1, ncol = 60))
names(merge_max_conRR) = c(1:60)
merge_max_conRR$second_category <- 'second_category'
merge_max_conRR$first_category <- 'first_category'




for(i in 1:n){         #对于每个一级目录(文件夹)
  b=list.files(dir[i]) #b是列出每个一级目录(文件夹)中每个xlsx文件的名称
  n_sub[i]=length(b)   #得到一级目录(文件夹)下xlsx的文件个数:n_sub
  

  for(j in 1:n_sub[i]){     #对于每个一级目录(文件夹)下的每个xlsx文件
    data<-read.delim(file=paste(dir[i],'/',b[j],sep=''),sep = "", header = F) #读取xlsx文件
    #conRRcum (mm)
    data_conRRcum = data[691:743,]
    data_conRRcum[,48:83]=trunc(data_conRRcum[,48:83]/1000000)
    data_conRRcum[2,64]=data_conRRcum[2,64]/10^234
    data_conRRcum = as.data.frame(lapply(data_conRRcum,as.numeric))
    data_conRRcum[2:53,]=data_conRRcum[2:53,]/10
    conRRcum_6hstep = data_conRRcum[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
    #conv precip (mm)
    conRR_6hstep = conRRcum_6hstep[,1:(ncol(conRRcum_6hstep)-1)]
    conRR_6hstep = cbind(0,conRR_6hstep)
    conRRcum_6hstep = as.data.frame(lapply(conRRcum_6hstep,as.numeric))
    conRR_6hstep = as.data.frame(lapply(conRR_6hstep,as.numeric))
    conRR_6hstep = conRRcum_6hstep-conRR_6hstep
    conRR_HRES = conRR_6hstep[2,1:40]
    conRR_ENS = conRR_6hstep[3:53,]
    conRR_ENS = as.data.frame(lapply(conRR_ENS,as.numeric))
    avg_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    sd_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    median_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q10_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q90_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q25_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q75_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    p0.05_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    p0.3_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    p10_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    min_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    max_var_df <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    for(k in 1:ncol(conRR_ENS)){ 
      #  p0 <- colSums(RR_ENS>0)/length(RR_ENS[,k])   
      avg_var_df[,k] <- mean(conRR_ENS[, k])
      sd_var_df[,k] <- sd(conRR_ENS[,k])
      median_var_df[,k] <- median(conRR_ENS[,k])
      q10_var_df[,k] = quantile(conRR_ENS[,k],probs=0.1)
      q25_var_df[,k] = quantile(conRR_ENS[,k],probs=0.25)
      q75_var_df[,k] = quantile(conRR_ENS[,k],probs=0.75)
      q90_var_df[,k] = quantile(conRR_ENS[,k],probs=0.9)
      p0.05_var_df[,k] = sum(conRR_ENS[,k] > 0.05)/length(conRR_ENS[,k])
      p0.3_var_df[,k] = sum(conRR_ENS[,k] > 0.3)/length(conRR_ENS[,k])
      p10_var_df[,k] = sum(conRR_ENS[,k] > 10)/length(conRR_ENS[,k])
      min_var_df[,k] <- min(conRR_ENS[, k])
      max_var_df[,k] <- max(conRR_ENS[,k])

    }    
    
    
    #new row for HRES
    #    new_HRES = RR_HRES
    names(conRR_HRES)<-c(1:40)
    conRR_HRES$second_category<-substr(b[j],14,27)      
    conRR_HRES$first_category<-first_category_name[i]   
    merge_HRES_conRR<-rbind(merge_HRES_conRR,conRR_HRES)
    #new row for avg
    #    new_avg = avg_var_df
    names(avg_var_df)<-c(1:60)
    avg_var_df$second_category<-substr(b[j],14,27)       
    avg_var_df$first_category<-first_category_name[i] 
    merge_avg_conRR<-rbind(merge_avg_conRR,avg_var_df)    
    #new row for sd
    #    new_sd = sd_var_df
    names(sd_var_df)<-c(1:60)
    sd_var_df$second_category<-substr(b[j],14,27)       
    sd_var_df$first_category<-first_category_name[i] 
    merge_sd_conRR<-rbind(merge_sd_conRR,sd_var_df) 
    #new row for median
    #    new_median = median_var_df
    names(median_var_df)<-c(1:60)
    median_var_df$second_category<-substr(b[j],14,27)       
    median_var_df$first_category<-first_category_name[i] 
    merge_median_conRR<-rbind(merge_median_conRR,median_var_df)
    #new row for q10
    #    new_q10 = q10_var_df
    names(q10_var_df)<-c(1:60)
    q10_var_df$second_category<-substr(b[j],14,27)       
    q10_var_df$first_category<-first_category_name[i] 
    merge_q10_conRR<-rbind(merge_q10_conRR,q10_var_df)
    #new row for q90
    #    new_q90 = q90_var_df
    names(q90_var_df)<-c(1:60)
    q90_var_df$second_category<-substr(b[j],14,27)       
    q90_var_df$first_category<-first_category_name[i] 
    merge_q90_conRR<-rbind(merge_q90_conRR,q90_var_df)
    #new row for q25
    names(q25_var_df)<-c(1:60)
    q25_var_df$second_category<-substr(b[j],14,27)       
    q25_var_df$first_category<-first_category_name[i] 
    merge_q25_conRR<-rbind(merge_q25_conRR,q25_var_df)
    #new row for q75
    names(q75_var_df)<-c(1:60)
    q75_var_df$second_category<-substr(b[j],14,27)       
    q75_var_df$first_category<-first_category_name[i] 
    merge_q75_conRR<-rbind(merge_q75_conRR,q75_var_df)
    
    #new row for p0
    #    new_p0 = p0_var_df
    names(p0.05_var_df)<-c(1:60)
    p0.05_var_df$second_category<-substr(b[j],14,27)       
    p0.05_var_df$first_category<-first_category_name[i] 
    merge_p0.05_conRR<-rbind(merge_p0.05_conRR,p0.05_var_df)    
    #new row for p1
    #    new_p1 = p1_var_df
    names(p0.3_var_df)<-c(1:60)
    p0.3_var_df$second_category<-substr(b[j],14,27)       
    p0.3_var_df$first_category<-first_category_name[i] 
    merge_p0.3_conRR<-rbind(merge_p0.3_conRR,p0.3_var_df)   
    
    #new row for p10
    #    new_p10 = p10_var_df
    names(p10_var_df)<-c(1:60)
    p10_var_df$second_category<-substr(b[j],14,27)       
    p10_var_df$first_category<-first_category_name[i] 
    merge_p10_conRR<-rbind(merge_p10_conRR,p10_var_df)  

    #new row for min
    names(min_var_df)<-c(1:60)
    min_var_df$second_category<-substr(b[j],14,27)       
    min_var_df$first_category<-first_category_name[i] 
    merge_min_conRR<-rbind(merge_min_conRR,min_var_df)    
    #new row for max
    names(max_var_df)<-c(1:60)
    max_var_df$second_category<-substr(b[j],14,27)       
    max_var_df$first_category<-first_category_name[i] 
    merge_max_conRR<-rbind(merge_max_conRR,max_var_df)    
  }
}    

Sys.time()
setwd("E:/ECMWFdata/pre-process")
write.table(merge_HRES_conRR, "conRR_6h_HRES_2018_12.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_avg_conRR, "conRR_6h_mean_2018_12.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_sd_conRR, "conRR_6h_sd_2018_12.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_median_conRR, "conRR_6h_median_2018_12.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q10_conRR, "conRR_6h_q10_2018_12.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q90_conRR, "conRR_6h_q90_2018_12.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q25_conRR, "conRR_6h_q25_2018_12.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q75_conRR, "conRR_6h_q75_2018_12.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_p0.05_conRR, "conRR_6h_p0.05_2018_12.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_p0.3_conRR, "conRR_6h_p0.3_2018_12.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_p10_conRR, "conRR_6h_p10_2018_12.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_min_conRR, "conRR_6h_min_2018_12.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_max_conRR, "conRR_6h_max_2018_12.txt", sep = "\t", col.names = T, row.names = F)


#######(For other predictors)#####################################
#create the df for avg
merge_avg_CC = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_avg_CC) = c(1:60)
merge_avg_CC$second_category <- 'second_category'
merge_avg_CC$first_category <- 'first_category'

merge_avg_cape = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_avg_cape) = c(1:60)
merge_avg_cape$second_category <- 'second_category'
merge_avg_cape$first_category <- 'first_category'

merge_avg_capeS = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_avg_capeS) = c(1:60)
merge_avg_capeS$second_category <- 'second_category'
merge_avg_capeS$first_category <- 'first_category'

merge_avg_eva = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_avg_eva) = c(1:60)
merge_avg_eva$second_category <- 'second_category'
merge_avg_eva$first_category <- 'first_category'

#create the df for sd
merge_sd_CC = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_sd_CC) = c(1:60)
merge_sd_CC$second_category <- 'second_category'
merge_sd_CC$first_category <- 'first_category'

merge_sd_cape = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_sd_cape) = c(1:60)
merge_sd_cape$second_category <- 'second_category'
merge_sd_cape$first_category <- 'first_category'

merge_sd_capeS = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_sd_capeS) = c(1:60)
merge_sd_capeS$second_category <- 'second_category'
merge_sd_capeS$first_category <- 'first_category'

merge_sd_eva = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_sd_eva) = c(1:60)
merge_sd_eva$second_category <- 'second_category'
merge_sd_eva$first_category <- 'first_category'

#create the df for median
merge_median_CC = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_median_CC) = c(1:60)
merge_median_CC$second_category <- 'second_category'
merge_median_CC$first_category <- 'first_category'

merge_median_cape = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_median_cape) = c(1:60)
merge_median_cape$second_category <- 'second_category'
merge_median_cape$first_category <- 'first_category'

merge_median_capeS = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_median_capeS) = c(1:60)
merge_median_capeS$second_category <- 'second_category'
merge_median_capeS$first_category <- 'first_category'

merge_median_eva = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_median_eva) = c(1:60)
merge_median_eva$second_category <- 'second_category'
merge_median_eva$first_category <- 'first_category'

#create the df for q10
merge_q10_CC = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q10_CC) = c(1:60)
merge_q10_CC$second_category <- 'second_category'
merge_q10_CC$first_category <- 'first_category'

merge_q10_cape = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q10_cape) = c(1:60)
merge_q10_cape$second_category <- 'second_category'
merge_q10_cape$first_category <- 'first_category'

merge_q10_capeS = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q10_capeS) = c(1:60)
merge_q10_capeS$second_category <- 'second_category'
merge_q10_capeS$first_category <- 'first_category'

merge_q10_eva = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q10_eva) = c(1:60)
merge_q10_eva$second_category <- 'second_category'
merge_q10_eva$first_category <- 'first_category'

#create the df for q25
merge_q25_CC = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q25_CC) = c(1:60)
merge_q25_CC$second_category <- 'second_category'
merge_q25_CC$first_category <- 'first_category'

merge_q25_cape = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q25_cape) = c(1:60)
merge_q25_cape$second_category <- 'second_category'
merge_q25_cape$first_category <- 'first_category'

merge_q25_capeS = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q25_capeS) = c(1:60)
merge_q25_capeS$second_category <- 'second_category'
merge_q25_capeS$first_category <- 'first_category'

merge_q25_eva = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q25_eva) = c(1:60)
merge_q25_eva$second_category <- 'second_category'
merge_q25_eva$first_category <- 'first_category'


#createe the df for q75
merge_q75_CC = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q75_CC) = c(1:60)
merge_q75_CC$second_category <- 'second_category'
merge_q75_CC$first_category <- 'first_category'

merge_q75_cape = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q75_cape) = c(1:60)
merge_q75_cape$second_category <- 'second_category'
merge_q75_cape$first_category <- 'first_category'

merge_q75_capeS = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q75_capeS) = c(1:60)
merge_q75_capeS$second_category <- 'second_category'
merge_q75_capeS$first_category <- 'first_category'

merge_q75_eva = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q75_eva) = c(1:60)
merge_q75_eva$second_category <- 'second_category'
merge_q75_eva$first_category <- 'first_category'

#create the df for q90
merge_q90_CC = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q90_CC) = c(1:60)
merge_q90_CC$second_category <- 'second_category'
merge_q90_CC$first_category <- 'first_category'

merge_q90_cape = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q90_cape) = c(1:60)
merge_q90_cape$second_category <- 'second_category'
merge_q90_cape$first_category <- 'first_category'

merge_q90_capeS = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q90_capeS) = c(1:60)
merge_q90_capeS$second_category <- 'second_category'
merge_q90_capeS$first_category <- 'first_category'

merge_q90_eva = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q90_eva) = c(1:60)
merge_q90_eva$second_category <- 'second_category'
merge_q90_eva$first_category <- 'first_category'

#create the df for min
merge_min_CC = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_min_CC) = c(1:60)
merge_min_CC$second_category <- 'second_category'
merge_min_CC$first_category <- 'first_category'

merge_min_cape = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_min_cape) = c(1:60)
merge_min_cape$second_category <- 'second_category'
merge_min_cape$first_category <- 'first_category'

merge_min_capeS = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_min_capeS) = c(1:60)
merge_min_capeS$second_category <- 'second_category'
merge_min_capeS$first_category <- 'first_category'

merge_min_eva = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_min_eva) = c(1:60)
merge_min_eva$second_category <- 'second_category'
merge_min_eva$first_category <- 'first_category'

#create the df for max
merge_max_CC = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_max_CC) = c(1:60)
merge_max_CC$second_category <- 'second_category'
merge_max_CC$first_category <- 'first_category'

merge_max_cape = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_max_cape) = c(1:60)
merge_max_cape$second_category <- 'second_category'
merge_max_cape$first_category <- 'first_category'

merge_max_capeS = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_max_capeS) = c(1:60)
merge_max_capeS$second_category <- 'second_category'
merge_max_capeS$first_category <- 'first_category'

merge_max_eva = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_max_eva) = c(1:60)
merge_max_eva$second_category <- 'second_category'
merge_max_eva$first_category <- 'first_category'



#####(For WSPD)#################################################################################################
#WSPD (m/s)
Sys.time()
setwd('E:/ECMWFdata/ECME_2021')
first_category_name = list.files("Mar")
dir = paste("./Mar/",first_category_name,sep="")
n = length(dir) 
n_sub<-rep(0,n)
n_sub<-as.data.frame(n_sub)
n_sub<-t(n_sub)
head(n_sub)
#create the df for HRES
merge_HRES_WSPD = as.data.frame(matrix(1:40,nrow = 1, ncol = 40))
names(merge_HRES_WSPD) = c(1:40)
merge_HRES_WSPD$second_category <- 'second_category'
merge_HRES_WSPD$first_category <- 'first_category'

merge_HRES_T2m = as.data.frame(matrix(1:40,nrow = 1, ncol = 40))
names(merge_HRES_T2m) = c(1:40)
merge_HRES_T2m$second_category <- 'second_category'
merge_HRES_T2m$first_category <- 'first_category'

merge_HRES_DPT = as.data.frame(matrix(1:40,nrow = 1, ncol = 40))
names(merge_HRES_DPT) = c(1:40)
merge_HRES_DPT$second_category <- 'second_category'
merge_HRES_DPT$first_category <- 'first_category'

#create the df for avg
merge_avg_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_avg_WSPD) = c(1:60)
merge_avg_WSPD$second_category <- 'second_category'
merge_avg_WSPD$first_category <- 'first_category'

merge_avg_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_avg_T2m) = c(1:60)
merge_avg_T2m$second_category <- 'second_category'
merge_avg_T2m$first_category <- 'first_category'

merge_avg_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_avg_DPT) = c(1:60)
merge_avg_DPT$second_category <- 'second_category'
merge_avg_DPT$first_category <- 'first_category'


#create the df for sd
merge_sd_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_sd_WSPD) = c(1:60)
merge_sd_WSPD$second_category <- 'second_category'
merge_sd_WSPD$first_category <- 'first_category'

merge_sd_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_sd_T2m) = c(1:60)
merge_sd_T2m$second_category <- 'second_category'
merge_sd_T2m$first_category <- 'first_category'

merge_sd_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_sd_DPT) = c(1:60)
merge_sd_DPT$second_category <- 'second_category'
merge_sd_DPT$first_category <- 'first_category'

#create the df for median
merge_median_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_median_WSPD) = c(1:60)
merge_median_WSPD$second_category <- 'second_category'
merge_median_WSPD$first_category <- 'first_category'

merge_median_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_median_T2m) = c(1:60)
merge_median_T2m$second_category <- 'second_category'
merge_median_T2m$first_category <- 'first_category'

merge_median_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_median_DPT) = c(1:60)
merge_median_DPT$second_category <- 'second_category'
merge_median_DPT$first_category <- 'first_category'

#create the df for q10
merge_q10_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q10_WSPD) = c(1:60)
merge_q10_WSPD$second_category <- 'second_category'
merge_q10_WSPD$first_category <- 'first_category'

merge_q10_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q10_T2m) = c(1:60)
merge_q10_T2m$second_category <- 'second_category'
merge_q10_T2m$first_category <- 'first_category'

merge_q10_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q10_DPT) = c(1:60)
merge_q10_DPT$second_category <- 'second_category'
merge_q10_DPT$first_category <- 'first_category'

#create the df for q25
merge_q25_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q25_WSPD) = c(1:60)
merge_q25_WSPD$second_category <- 'second_category'
merge_q25_WSPD$first_category <- 'first_category'

merge_q25_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q25_T2m) = c(1:60)
merge_q25_T2m$second_category <- 'second_category'
merge_q25_T2m$first_category <- 'first_category'

merge_q25_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q25_DPT) = c(1:60)
merge_q25_DPT$second_category <- 'second_category'
merge_q25_DPT$first_category <- 'first_category'

#createe the df for q75
merge_q75_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q75_WSPD) = c(1:60)
merge_q75_WSPD$second_category <- 'second_category'
merge_q75_WSPD$first_category <- 'first_category'

merge_q75_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q75_T2m) = c(1:60)
merge_q75_T2m$second_category <- 'second_category'
merge_q75_T2m$first_category <- 'first_category'

merge_q75_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q75_DPT) = c(1:60)
merge_q75_DPT$second_category <- 'second_category'
merge_q75_DPT$first_category <- 'first_category'


#create the df for q90
merge_q90_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q90_WSPD) = c(1:60)
merge_q90_WSPD$second_category <- 'second_category'
merge_q90_WSPD$first_category <- 'first_category'

merge_q90_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q90_T2m) = c(1:60)
merge_q90_T2m$second_category <- 'second_category'
merge_q90_T2m$first_category <- 'first_category'

merge_q90_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q90_DPT) = c(1:60)
merge_q90_DPT$second_category <- 'second_category'
merge_q90_DPT$first_category <- 'first_category'


#create the df for min
merge_min_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_min_WSPD) = c(1:60)
merge_min_WSPD$second_category <- 'second_category'
merge_min_WSPD$first_category <- 'first_category'

merge_min_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_min_T2m) = c(1:60)
merge_min_T2m$second_category <- 'second_category'
merge_min_T2m$first_category <- 'first_category'

merge_min_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_min_DPT) = c(1:60)
merge_min_DPT$second_category <- 'second_category'
merge_min_DPT$first_category <- 'first_category'


#create the df for max
merge_max_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_max_WSPD) = c(1:60)
merge_max_WSPD$second_category <- 'second_category'
merge_max_WSPD$first_category <- 'first_category'

merge_max_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_max_T2m) = c(1:60)
merge_max_T2m$second_category <- 'second_category'
merge_max_T2m$first_category <- 'first_category'

merge_max_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_max_DPT) = c(1:60)
merge_max_DPT$second_category <- 'second_category'
merge_max_DPT$first_category <- 'first_category'






for(i in 1:n){         #对于每个一级目录(文件夹)
  b=list.files(dir[i]) #b是列出每个一级目录(文件夹)中每个xlsx文件的名称
  n_sub[i]=length(b)   #得到一级目录(文件夹)下xlsx的文件个数:n_sub
  
  for(j in 1:n_sub[i]){     #对于每个一级目录(文件夹)下的每个xlsx文件
    data<-read.delim(file=paste(dir[i],'/',b[j],sep=''),sep = "", header = F) #读取xlsx文件
    data_WSPD = data[161:213,]
    data_WSPD[,48:83]=trunc(data_WSPD[,48:83]/1000000)
    data_WSPD[2,64]=data_WSPD[2,64]/10^234
    data_WSPD = as.data.frame(lapply(data_WSPD,as.numeric))
    data_WSPD[2:53,]=data_WSPD[2:53,]/10
    WSPD_6hstep = data_WSPD[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
    #T2m (degrees Celcius)
    data_T2m = data[267:319,]
    data_T2m[,48:83]=trunc(data_T2m[,48:83]/1000000)
    data_T2m[2,64]=data_T2m[2,64]/10^234
    data_T2m = as.data.frame(lapply(data_T2m,as.numeric))
    data_T2m[2:53,]=data_T2m[2:53,]/10
    T2m_6hstep = data_T2m[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
    #DPT (degrees Celcius)
    data_DPT = data[320:372,]
    data_DPT[,48:83]=trunc(data_DPT[,48:83]/1000000)
    data_DPT[2,64]=data_DPT[2,64]/10^234
    data_DPT = as.data.frame(lapply(data_DPT,as.numeric))
    data_DPT[2:53,]=data_DPT[2:53,]/10
    DPT_6hstep = data_DPT[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
   
    WSPD_HRES = WSPD_6hstep[2,1:40]
    DPT_HRES = DPT_6hstep[2,1:40]
    T2m_HRES = T2m_6hstep[2,1:40]
    
    WSPD_ENS = WSPD_6hstep[3:53,]
    WSPD_ENS = as.data.frame(lapply(WSPD_ENS,as.numeric))
    
    T2m_ENS = T2m_6hstep[3:53,]
    T2m_ENS = as.data.frame(lapply(T2m_ENS,as.numeric))
    
    DPT_ENS = DPT_6hstep[3:53,]
    DPT_ENS = as.data.frame(lapply(DPT_ENS,as.numeric))
    
    avg_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    avg_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    avg_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    sd_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    sd_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    sd_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    
    median_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    median_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    median_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    
    q10_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q10_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q10_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    q25_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q25_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q25_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    q75_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q75_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q75_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    
    q90_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q90_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q90_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
  
    min_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    min_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    min_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    max_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    max_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    max_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    
    
    for(k in 1:60){ 
      avg_var_df_WSPD[,k] <- mean(WSPD_ENS[, k])
      avg_var_df_T2m[,k] <- mean(T2m_ENS[, k])
      avg_var_df_DPT[,k] <- mean(DPT_ENS[, k])
      
      sd_var_df_WSPD[,k] <- sd(WSPD_ENS[, k])
      sd_var_df_T2m[,k] <- sd(T2m_ENS[, k])
      sd_var_df_DPT[,k] <- sd(DPT_ENS[, k])
      
      median_var_df_WSPD[,k] <- median(WSPD_ENS[, k])
      median_var_df_T2m[,k] <- median(T2m_ENS[, k])
      median_var_df_DPT[,k] <- median(DPT_ENS[, k])
      
      q10_var_df_WSPD[,k] <- quantile(WSPD_ENS[,k],probs=0.1)
      q10_var_df_T2m[,k] <- quantile(T2m_ENS[,k],probs=0.1)
      q10_var_df_DPT[,k] <- quantile(DPT_ENS[,k],probs=0.1)
      
      q25_var_df_WSPD[,k] <- quantile(WSPD_ENS[,k],probs=0.25)
      q25_var_df_T2m[,k] <- quantile(T2m_ENS[,k],probs=0.25)
      q25_var_df_DPT[,k] <- quantile(DPT_ENS[,k],probs=0.25)
      
      q75_var_df_WSPD[,k] <- quantile(WSPD_ENS[,k],probs=0.75)
      q75_var_df_T2m[,k] <- quantile(T2m_ENS[,k],probs=0.75)
      q75_var_df_DPT[,k] <- quantile(DPT_ENS[,k],probs=0.75)
      
      q90_var_df_WSPD[,k] <- quantile(WSPD_ENS[,k],probs=0.9)
      q90_var_df_T2m[,k] <- quantile(T2m_ENS[,k],probs=0.9)
      q90_var_df_DPT[,k] <- quantile(DPT_ENS[,k],probs=0.9)
      
      
      min_var_df_WSPD[,k] <- min(WSPD_ENS[, k])
      min_var_df_T2m[,k] <- min(T2m_ENS[, k])
      min_var_df_DPT[,k] <- min(DPT_ENS[, k])
      
      max_var_df_WSPD[,k] <- max(WSPD_ENS[, k])
      max_var_df_T2m[,k] <- max(T2m_ENS[, k])
      max_var_df_DPT[,k] <- max(DPT_ENS[, k])
      
    }    
    #new row for HRES
    names(WSPD_HRES)<-c(1:40)
    WSPD_HRES$second_category<-substr(b[j],14,27)      
    WSPD_HRES$first_category<-first_category_name[i]   
    merge_HRES_WSPD<-rbind(merge_HRES_WSPD,WSPD_HRES)
    
    names(T2m_HRES)<-c(1:40)
    T2m_HRES$second_category<-substr(b[j],14,27)      
    T2m_HRES$first_category<-first_category_name[i]   
    merge_HRES_T2m<-rbind(merge_HRES_T2m,T2m_HRES)
    
    names(DPT_HRES)<-c(1:40)
    DPT_HRES$second_category<-substr(b[j],14,27)      
    DPT_HRES$first_category<-first_category_name[i]   
    merge_HRES_DPT<-rbind(merge_HRES_DPT,DPT_HRES)
    
   #new row for avg
    names(avg_var_df_WSPD)<-c(1:60)
    avg_var_df_WSPD$second_category<-substr(b[j],14,27)       
    avg_var_df_WSPD$first_category<-first_category_name[i] 
    merge_avg_WSPD<-rbind(merge_avg_WSPD,avg_var_df_WSPD)
    
    names(avg_var_df_T2m)<-c(1:60)
    avg_var_df_T2m$second_category<-substr(b[j],14,27)       
    avg_var_df_T2m$first_category<-first_category_name[i] 
    merge_avg_T2m<-rbind(merge_avg_T2m,avg_var_df_T2m)
    
    names(avg_var_df_DPT)<-c(1:60)
    avg_var_df_DPT$second_category<-substr(b[j],14,27)       
    avg_var_df_DPT$first_category<-first_category_name[i] 
    merge_avg_DPT<-rbind(merge_avg_DPT,avg_var_df_DPT)
    
    #new row for sd
    names(sd_var_df_WSPD)<-c(1:60)
    sd_var_df_WSPD$second_category<-substr(b[j],14,27)       
    sd_var_df_WSPD$first_category<-first_category_name[i] 
    merge_sd_WSPD<-rbind(merge_sd_WSPD,sd_var_df_WSPD)
    
    names(sd_var_df_T2m)<-c(1:60)
    sd_var_df_T2m$second_category<-substr(b[j],14,27)       
    sd_var_df_T2m$first_category<-first_category_name[i] 
    merge_sd_T2m<-rbind(merge_sd_T2m,sd_var_df_T2m)
    
    names(sd_var_df_DPT)<-c(1:60)
    sd_var_df_DPT$second_category<-substr(b[j],14,27)       
    sd_var_df_DPT$first_category<-first_category_name[i] 
    merge_sd_DPT<-rbind(merge_sd_DPT,sd_var_df_DPT)
    
    
    #new row for median
    names(median_var_df_WSPD)<-c(1:60)
    median_var_df_WSPD$second_category<-substr(b[j],14,27)       
    median_var_df_WSPD$first_category<-first_category_name[i] 
    merge_median_WSPD<-rbind(merge_median_WSPD,median_var_df_WSPD)
    
    names(median_var_df_T2m)<-c(1:60)
    median_var_df_T2m$second_category<-substr(b[j],14,27)       
    median_var_df_T2m$first_category<-first_category_name[i] 
    merge_median_T2m<-rbind(merge_median_T2m,median_var_df_T2m)
    
    names(median_var_df_DPT)<-c(1:60)
    median_var_df_DPT$second_category<-substr(b[j],14,27)       
    median_var_df_DPT$first_category<-first_category_name[i] 
    merge_median_DPT<-rbind(merge_median_DPT,median_var_df_DPT)
    
    #new row for q10
    names(q10_var_df_WSPD)<-c(1:60)
    q10_var_df_WSPD$second_category<-substr(b[j],14,27)       
    q10_var_df_WSPD$first_category<-first_category_name[i] 
    merge_q10_WSPD<-rbind(merge_q10_WSPD,q10_var_df_WSPD)
    
    names(q10_var_df_T2m)<-c(1:60)
    q10_var_df_T2m$second_category<-substr(b[j],14,27)       
    q10_var_df_T2m$first_category<-first_category_name[i] 
    merge_q10_T2m<-rbind(merge_q10_T2m,q10_var_df_T2m)
    
    names(q10_var_df_DPT)<-c(1:60)
    q10_var_df_DPT$second_category<-substr(b[j],14,27)       
    q10_var_df_DPT$first_category<-first_category_name[i] 
    merge_q10_DPT<-rbind(merge_q10_DPT,q10_var_df_DPT)
    
    #new row for q25
    names(q25_var_df_WSPD)<-c(1:60)
    q25_var_df_WSPD$second_category<-substr(b[j],14,27)       
    q25_var_df_WSPD$first_category<-first_category_name[i] 
    merge_q25_WSPD<-rbind(merge_q25_WSPD,q25_var_df_WSPD)
    
    names(q25_var_df_T2m)<-c(1:60)
    q25_var_df_T2m$second_category<-substr(b[j],14,27)       
    q25_var_df_T2m$first_category<-first_category_name[i] 
    merge_q25_T2m<-rbind(merge_q25_T2m,q25_var_df_T2m)
    
    names(q25_var_df_DPT)<-c(1:60)
    q25_var_df_DPT$second_category<-substr(b[j],14,27)       
    q25_var_df_DPT$first_category<-first_category_name[i] 
    merge_q25_DPT<-rbind(merge_q25_DPT,q25_var_df_DPT)
    
    #new row for q75
    names(q75_var_df_WSPD)<-c(1:60)
    q75_var_df_WSPD$second_category<-substr(b[j],14,27)       
    q75_var_df_WSPD$first_category<-first_category_name[i] 
    merge_q75_WSPD<-rbind(merge_q75_WSPD,q75_var_df_WSPD)
    
    names(q75_var_df_T2m)<-c(1:60)
    q75_var_df_T2m$second_category<-substr(b[j],14,27)       
    q75_var_df_T2m$first_category<-first_category_name[i] 
    merge_q75_T2m<-rbind(merge_q75_T2m,q75_var_df_T2m)
    
    names(q75_var_df_DPT)<-c(1:60)
    q75_var_df_DPT$second_category<-substr(b[j],14,27)       
    q75_var_df_DPT$first_category<-first_category_name[i] 
    merge_q75_DPT<-rbind(merge_q75_DPT,q75_var_df_DPT)
    
    
    
    #new row for q90
    names(q90_var_df_WSPD)<-c(1:60)
    q90_var_df_WSPD$second_category<-substr(b[j],14,27)       
    q90_var_df_WSPD$first_category<-first_category_name[i] 
    merge_q90_WSPD<-rbind(merge_q90_WSPD,q90_var_df_WSPD)
    
    names(q90_var_df_T2m)<-c(1:60)
    q90_var_df_T2m$second_category<-substr(b[j],14,27)       
    q90_var_df_T2m$first_category<-first_category_name[i] 
    merge_q90_T2m<-rbind(merge_q90_T2m,q90_var_df_T2m)
    
    names(q90_var_df_DPT)<-c(1:60)
    q90_var_df_DPT$second_category<-substr(b[j],14,27)       
    q90_var_df_DPT$first_category<-first_category_name[i] 
    merge_q90_DPT<-rbind(merge_q90_DPT,q90_var_df_DPT)
    
    #new row for min
    names(min_var_df_WSPD)<-c(1:60)
    min_var_df_WSPD$second_category<-substr(b[j],14,27)       
    min_var_df_WSPD$first_category<-first_category_name[i] 
    merge_min_WSPD<-rbind(merge_min_WSPD,min_var_df_WSPD)
    
    names(min_var_df_T2m)<-c(1:60)
    min_var_df_T2m$second_category<-substr(b[j],14,27)       
    min_var_df_T2m$first_category<-first_category_name[i] 
    merge_min_T2m<-rbind(merge_min_T2m,min_var_df_T2m)
    
    names(min_var_df_DPT)<-c(1:60)
    min_var_df_DPT$second_category<-substr(b[j],14,27)       
    min_var_df_DPT$first_category<-first_category_name[i] 
   merge_min_DPT<-rbind(merge_min_DPT,min_var_df_DPT)
    
    #new row for max
    names(max_var_df_WSPD)<-c(1:60)
    max_var_df_WSPD$second_category<-substr(b[j],14,27)       
    max_var_df_WSPD$first_category<-first_category_name[i] 
    merge_max_WSPD<-rbind(merge_max_WSPD,max_var_df_WSPD)
    
    names(max_var_df_T2m)<-c(1:60)
    max_var_df_T2m$second_category<-substr(b[j],14,27)       
    max_var_df_T2m$first_category<-first_category_name[i] 
    merge_max_T2m<-rbind(merge_max_T2m,max_var_df_T2m)
    
    names(max_var_df_DPT)<-c(1:60)
    max_var_df_DPT$second_category<-substr(b[j],14,27)       
    max_var_df_DPT$first_category<-first_category_name[i] 
    merge_max_DPT<-rbind(merge_max_DPT,max_var_df_DPT)
    
  }
}   
Sys.time()
setwd("E:/ECMWFdata/pre-process")
write.table(merge_HRES_WSPD, "WSPD_6h_HRES_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_avg_WSPD, "WSPD_6h_mean_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_sd_WSPD, "WSPD_6h_sd_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_median_WSPD, "WSPD_6h_median_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q10_WSPD, "WSPD_6h_q10_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q25_WSPD, "WSPD_6h_q25_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q75_WSPD, "WSPD_6h_q75_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q90_WSPD, "WSPD_6h_q90_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_min_WSPD, "WSPD_6h_min_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_max_WSPD, "WSPD_6h_max_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)

write.table(merge_HRES_T2m, "T2m_6h_HRES_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_avg_T2m, "T2m_6h_mean_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_sd_T2m, "T2m_6h_sd_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_median_T2m, "T2m_6h_median_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q10_T2m, "T2m_6h_q10_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q25_T2m, "T2m_6h_q25_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q75_T2m, "T2m_6h_q75_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q90_T2m, "T2m_6h_q90_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_min_T2m, "T2m_6h_min_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_max_T2m, "T2m_6h_max_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)

write.table(merge_HRES_DPT, "DPT_6h_HRES_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_avg_DPT, "DPT_6h_mean_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_sd_DPT, "DPT_6h_sd_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_median_DPT, "DPT_6h_median_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q10_DPT, "DPT_6h_q10_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q25_DPT, "DPT_6h_q25_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q75_DPT, "DPT_6h_q75_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q90_DPT, "DPT_6h_q90_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_min_DPT, "DPT_6h_min_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_max_DPT, "DPT_6h_max_2021_Mar.txt", sep = "\t", col.names = T, row.names = F)



rm(list = ls())




#WSPD (m/s)
Sys.time()
setwd('E:/ECMWFdata/ECME_2021')
first_category_name = list.files("Apr")
dir = paste("./Apr/",first_category_name,sep="")
n = length(dir) 
n_sub<-rep(0,n)
n_sub<-as.data.frame(n_sub)
n_sub<-t(n_sub)
head(n_sub)
#create the df for HRES
merge_HRES_WSPD = as.data.frame(matrix(1:40,nrow = 1, ncol = 40))
names(merge_HRES_WSPD) = c(1:40)
merge_HRES_WSPD$second_category <- 'second_category'
merge_HRES_WSPD$first_category <- 'first_category'

merge_HRES_T2m = as.data.frame(matrix(1:40,nrow = 1, ncol = 40))
names(merge_HRES_T2m) = c(1:40)
merge_HRES_T2m$second_category <- 'second_category'
merge_HRES_T2m$first_category <- 'first_category'

merge_HRES_DPT = as.data.frame(matrix(1:40,nrow = 1, ncol = 40))
names(merge_HRES_DPT) = c(1:40)
merge_HRES_DPT$second_category <- 'second_category'
merge_HRES_DPT$first_category <- 'first_category'

#create the df for avg
merge_avg_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_avg_WSPD) = c(1:60)
merge_avg_WSPD$second_category <- 'second_category'
merge_avg_WSPD$first_category <- 'first_category'

merge_avg_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_avg_T2m) = c(1:60)
merge_avg_T2m$second_category <- 'second_category'
merge_avg_T2m$first_category <- 'first_category'

merge_avg_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_avg_DPT) = c(1:60)
merge_avg_DPT$second_category <- 'second_category'
merge_avg_DPT$first_category <- 'first_category'


#create the df for sd
merge_sd_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_sd_WSPD) = c(1:60)
merge_sd_WSPD$second_category <- 'second_category'
merge_sd_WSPD$first_category <- 'first_category'

merge_sd_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_sd_T2m) = c(1:60)
merge_sd_T2m$second_category <- 'second_category'
merge_sd_T2m$first_category <- 'first_category'

merge_sd_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_sd_DPT) = c(1:60)
merge_sd_DPT$second_category <- 'second_category'
merge_sd_DPT$first_category <- 'first_category'

#create the df for median
merge_median_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_median_WSPD) = c(1:60)
merge_median_WSPD$second_category <- 'second_category'
merge_median_WSPD$first_category <- 'first_category'

merge_median_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_median_T2m) = c(1:60)
merge_median_T2m$second_category <- 'second_category'
merge_median_T2m$first_category <- 'first_category'

merge_median_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_median_DPT) = c(1:60)
merge_median_DPT$second_category <- 'second_category'
merge_median_DPT$first_category <- 'first_category'

#create the df for q10
merge_q10_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q10_WSPD) = c(1:60)
merge_q10_WSPD$second_category <- 'second_category'
merge_q10_WSPD$first_category <- 'first_category'

merge_q10_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q10_T2m) = c(1:60)
merge_q10_T2m$second_category <- 'second_category'
merge_q10_T2m$first_category <- 'first_category'

merge_q10_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q10_DPT) = c(1:60)
merge_q10_DPT$second_category <- 'second_category'
merge_q10_DPT$first_category <- 'first_category'

#create the df for q25
merge_q25_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q25_WSPD) = c(1:60)
merge_q25_WSPD$second_category <- 'second_category'
merge_q25_WSPD$first_category <- 'first_category'

merge_q25_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q25_T2m) = c(1:60)
merge_q25_T2m$second_category <- 'second_category'
merge_q25_T2m$first_category <- 'first_category'

merge_q25_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q25_DPT) = c(1:60)
merge_q25_DPT$second_category <- 'second_category'
merge_q25_DPT$first_category <- 'first_category'

#createe the df for q75
merge_q75_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q75_WSPD) = c(1:60)
merge_q75_WSPD$second_category <- 'second_category'
merge_q75_WSPD$first_category <- 'first_category'

merge_q75_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q75_T2m) = c(1:60)
merge_q75_T2m$second_category <- 'second_category'
merge_q75_T2m$first_category <- 'first_category'

merge_q75_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q75_DPT) = c(1:60)
merge_q75_DPT$second_category <- 'second_category'
merge_q75_DPT$first_category <- 'first_category'


#create the df for q90
merge_q90_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q90_WSPD) = c(1:60)
merge_q90_WSPD$second_category <- 'second_category'
merge_q90_WSPD$first_category <- 'first_category'

merge_q90_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q90_T2m) = c(1:60)
merge_q90_T2m$second_category <- 'second_category'
merge_q90_T2m$first_category <- 'first_category'

merge_q90_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q90_DPT) = c(1:60)
merge_q90_DPT$second_category <- 'second_category'
merge_q90_DPT$first_category <- 'first_category'


#create the df for min
merge_min_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_min_WSPD) = c(1:60)
merge_min_WSPD$second_category <- 'second_category'
merge_min_WSPD$first_category <- 'first_category'

merge_min_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_min_T2m) = c(1:60)
merge_min_T2m$second_category <- 'second_category'
merge_min_T2m$first_category <- 'first_category'

merge_min_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_min_DPT) = c(1:60)
merge_min_DPT$second_category <- 'second_category'
merge_min_DPT$first_category <- 'first_category'


#create the df for max
merge_max_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_max_WSPD) = c(1:60)
merge_max_WSPD$second_category <- 'second_category'
merge_max_WSPD$first_category <- 'first_category'

merge_max_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_max_T2m) = c(1:60)
merge_max_T2m$second_category <- 'second_category'
merge_max_T2m$first_category <- 'first_category'

merge_max_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_max_DPT) = c(1:60)
merge_max_DPT$second_category <- 'second_category'
merge_max_DPT$first_category <- 'first_category'






for(i in 1:n){         #对于每个一级目录(文件夹)
  b=list.files(dir[i]) #b是列出每个一级目录(文件夹)中每个xlsx文件的名称
  n_sub[i]=length(b)   #得到一级目录(文件夹)下xlsx的文件个数:n_sub
  
  for(j in 1:n_sub[i]){     #对于每个一级目录(文件夹)下的每个xlsx文件
    data<-read.delim(file=paste(dir[i],'/',b[j],sep=''),sep = "", header = F) #读取xlsx文件
    data_WSPD = data[161:213,]
    data_WSPD[,48:83]=trunc(data_WSPD[,48:83]/1000000)
    data_WSPD[2,64]=data_WSPD[2,64]/10^234
    data_WSPD = as.data.frame(lapply(data_WSPD,as.numeric))
    data_WSPD[2:53,]=data_WSPD[2:53,]/10
    WSPD_6hstep = data_WSPD[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
    #T2m (degrees Celcius)
    data_T2m = data[267:319,]
    data_T2m[,48:83]=trunc(data_T2m[,48:83]/1000000)
    data_T2m[2,64]=data_T2m[2,64]/10^234
    data_T2m = as.data.frame(lapply(data_T2m,as.numeric))
    data_T2m[2:53,]=data_T2m[2:53,]/10
    T2m_6hstep = data_T2m[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
    #DPT (degrees Celcius)
    data_DPT = data[320:372,]
    data_DPT[,48:83]=trunc(data_DPT[,48:83]/1000000)
    data_DPT[2,64]=data_DPT[2,64]/10^234
    data_DPT = as.data.frame(lapply(data_DPT,as.numeric))
    data_DPT[2:53,]=data_DPT[2:53,]/10
    DPT_6hstep = data_DPT[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
    
    WSPD_HRES = WSPD_6hstep[2,1:40]
    DPT_HRES = DPT_6hstep[2,1:40]
    T2m_HRES = T2m_6hstep[2,1:40]
    
    WSPD_ENS = WSPD_6hstep[3:53,]
    WSPD_ENS = as.data.frame(lapply(WSPD_ENS,as.numeric))
    
    T2m_ENS = T2m_6hstep[3:53,]
    T2m_ENS = as.data.frame(lapply(T2m_ENS,as.numeric))
    
    DPT_ENS = DPT_6hstep[3:53,]
    DPT_ENS = as.data.frame(lapply(DPT_ENS,as.numeric))
    
    avg_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    avg_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    avg_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    sd_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    sd_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    sd_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    
    median_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    median_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    median_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    
    q10_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q10_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q10_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    q25_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q25_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q25_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    q75_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q75_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q75_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    
    q90_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q90_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q90_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    min_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    min_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    min_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    max_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    max_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    max_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    
    
    for(k in 1:60){ 
      avg_var_df_WSPD[,k] <- mean(WSPD_ENS[, k])
      avg_var_df_T2m[,k] <- mean(T2m_ENS[, k])
      avg_var_df_DPT[,k] <- mean(DPT_ENS[, k])
      
      sd_var_df_WSPD[,k] <- sd(WSPD_ENS[, k])
      sd_var_df_T2m[,k] <- sd(T2m_ENS[, k])
      sd_var_df_DPT[,k] <- sd(DPT_ENS[, k])
      
      median_var_df_WSPD[,k] <- median(WSPD_ENS[, k])
      median_var_df_T2m[,k] <- median(T2m_ENS[, k])
      median_var_df_DPT[,k] <- median(DPT_ENS[, k])
      
      q10_var_df_WSPD[,k] <- quantile(WSPD_ENS[,k],probs=0.1)
      q10_var_df_T2m[,k] <- quantile(T2m_ENS[,k],probs=0.1)
      q10_var_df_DPT[,k] <- quantile(DPT_ENS[,k],probs=0.1)
      
      q25_var_df_WSPD[,k] <- quantile(WSPD_ENS[,k],probs=0.25)
      q25_var_df_T2m[,k] <- quantile(T2m_ENS[,k],probs=0.25)
      q25_var_df_DPT[,k] <- quantile(DPT_ENS[,k],probs=0.25)
      
      q75_var_df_WSPD[,k] <- quantile(WSPD_ENS[,k],probs=0.75)
      q75_var_df_T2m[,k] <- quantile(T2m_ENS[,k],probs=0.75)
      q75_var_df_DPT[,k] <- quantile(DPT_ENS[,k],probs=0.75)
      
      q90_var_df_WSPD[,k] <- quantile(WSPD_ENS[,k],probs=0.9)
      q90_var_df_T2m[,k] <- quantile(T2m_ENS[,k],probs=0.9)
      q90_var_df_DPT[,k] <- quantile(DPT_ENS[,k],probs=0.9)
      
      
      min_var_df_WSPD[,k] <- min(WSPD_ENS[, k])
      min_var_df_T2m[,k] <- min(T2m_ENS[, k])
      min_var_df_DPT[,k] <- min(DPT_ENS[, k])
      
      max_var_df_WSPD[,k] <- max(WSPD_ENS[, k])
      max_var_df_T2m[,k] <- max(T2m_ENS[, k])
      max_var_df_DPT[,k] <- max(DPT_ENS[, k])
      
    }    
    #new row for HRES
    names(WSPD_HRES)<-c(1:40)
    WSPD_HRES$second_category<-substr(b[j],14,27)      
    WSPD_HRES$first_category<-first_category_name[i]   
    merge_HRES_WSPD<-rbind(merge_HRES_WSPD,WSPD_HRES)
    
    names(T2m_HRES)<-c(1:40)
    T2m_HRES$second_category<-substr(b[j],14,27)      
    T2m_HRES$first_category<-first_category_name[i]   
    merge_HRES_T2m<-rbind(merge_HRES_T2m,T2m_HRES)
    
    names(DPT_HRES)<-c(1:40)
    DPT_HRES$second_category<-substr(b[j],14,27)      
    DPT_HRES$first_category<-first_category_name[i]   
    merge_HRES_DPT<-rbind(merge_HRES_DPT,DPT_HRES)
    
    #new row for avg
    names(avg_var_df_WSPD)<-c(1:60)
    avg_var_df_WSPD$second_category<-substr(b[j],14,27)       
    avg_var_df_WSPD$first_category<-first_category_name[i] 
    merge_avg_WSPD<-rbind(merge_avg_WSPD,avg_var_df_WSPD)
    
    names(avg_var_df_T2m)<-c(1:60)
    avg_var_df_T2m$second_category<-substr(b[j],14,27)       
    avg_var_df_T2m$first_category<-first_category_name[i] 
    merge_avg_T2m<-rbind(merge_avg_T2m,avg_var_df_T2m)
    
    names(avg_var_df_DPT)<-c(1:60)
    avg_var_df_DPT$second_category<-substr(b[j],14,27)       
    avg_var_df_DPT$first_category<-first_category_name[i] 
    merge_avg_DPT<-rbind(merge_avg_DPT,avg_var_df_DPT)
    
    #new row for sd
    names(sd_var_df_WSPD)<-c(1:60)
    sd_var_df_WSPD$second_category<-substr(b[j],14,27)       
    sd_var_df_WSPD$first_category<-first_category_name[i] 
    merge_sd_WSPD<-rbind(merge_sd_WSPD,sd_var_df_WSPD)
    
    names(sd_var_df_T2m)<-c(1:60)
    sd_var_df_T2m$second_category<-substr(b[j],14,27)       
    sd_var_df_T2m$first_category<-first_category_name[i] 
    merge_sd_T2m<-rbind(merge_sd_T2m,sd_var_df_T2m)
    
    names(sd_var_df_DPT)<-c(1:60)
    sd_var_df_DPT$second_category<-substr(b[j],14,27)       
    sd_var_df_DPT$first_category<-first_category_name[i] 
    merge_sd_DPT<-rbind(merge_sd_DPT,sd_var_df_DPT)
    
    
    #new row for median
    names(median_var_df_WSPD)<-c(1:60)
    median_var_df_WSPD$second_category<-substr(b[j],14,27)       
    median_var_df_WSPD$first_category<-first_category_name[i] 
    merge_median_WSPD<-rbind(merge_median_WSPD,median_var_df_WSPD)
    
    names(median_var_df_T2m)<-c(1:60)
    median_var_df_T2m$second_category<-substr(b[j],14,27)       
    median_var_df_T2m$first_category<-first_category_name[i] 
    merge_median_T2m<-rbind(merge_median_T2m,median_var_df_T2m)
    
    names(median_var_df_DPT)<-c(1:60)
    median_var_df_DPT$second_category<-substr(b[j],14,27)       
    median_var_df_DPT$first_category<-first_category_name[i] 
    merge_median_DPT<-rbind(merge_median_DPT,median_var_df_DPT)
    
    #new row for q10
    names(q10_var_df_WSPD)<-c(1:60)
    q10_var_df_WSPD$second_category<-substr(b[j],14,27)       
    q10_var_df_WSPD$first_category<-first_category_name[i] 
    merge_q10_WSPD<-rbind(merge_q10_WSPD,q10_var_df_WSPD)
    
    names(q10_var_df_T2m)<-c(1:60)
    q10_var_df_T2m$second_category<-substr(b[j],14,27)       
    q10_var_df_T2m$first_category<-first_category_name[i] 
    merge_q10_T2m<-rbind(merge_q10_T2m,q10_var_df_T2m)
    
    names(q10_var_df_DPT)<-c(1:60)
    q10_var_df_DPT$second_category<-substr(b[j],14,27)       
    q10_var_df_DPT$first_category<-first_category_name[i] 
    merge_q10_DPT<-rbind(merge_q10_DPT,q10_var_df_DPT)
    
    #new row for q25
    names(q25_var_df_WSPD)<-c(1:60)
    q25_var_df_WSPD$second_category<-substr(b[j],14,27)       
    q25_var_df_WSPD$first_category<-first_category_name[i] 
    merge_q25_WSPD<-rbind(merge_q25_WSPD,q25_var_df_WSPD)
    
    names(q25_var_df_T2m)<-c(1:60)
    q25_var_df_T2m$second_category<-substr(b[j],14,27)       
    q25_var_df_T2m$first_category<-first_category_name[i] 
    merge_q25_T2m<-rbind(merge_q25_T2m,q25_var_df_T2m)
    
    names(q25_var_df_DPT)<-c(1:60)
    q25_var_df_DPT$second_category<-substr(b[j],14,27)       
    q25_var_df_DPT$first_category<-first_category_name[i] 
    merge_q25_DPT<-rbind(merge_q25_DPT,q25_var_df_DPT)
    
    #new row for q75
    names(q75_var_df_WSPD)<-c(1:60)
    q75_var_df_WSPD$second_category<-substr(b[j],14,27)       
    q75_var_df_WSPD$first_category<-first_category_name[i] 
    merge_q75_WSPD<-rbind(merge_q75_WSPD,q75_var_df_WSPD)
    
    names(q75_var_df_T2m)<-c(1:60)
    q75_var_df_T2m$second_category<-substr(b[j],14,27)       
    q75_var_df_T2m$first_category<-first_category_name[i] 
    merge_q75_T2m<-rbind(merge_q75_T2m,q75_var_df_T2m)
    
    names(q75_var_df_DPT)<-c(1:60)
    q75_var_df_DPT$second_category<-substr(b[j],14,27)       
    q75_var_df_DPT$first_category<-first_category_name[i] 
    merge_q75_DPT<-rbind(merge_q75_DPT,q75_var_df_DPT)
    
    
    
    #new row for q90
    names(q90_var_df_WSPD)<-c(1:60)
    q90_var_df_WSPD$second_category<-substr(b[j],14,27)       
    q90_var_df_WSPD$first_category<-first_category_name[i] 
    merge_q90_WSPD<-rbind(merge_q90_WSPD,q90_var_df_WSPD)
    
    names(q90_var_df_T2m)<-c(1:60)
    q90_var_df_T2m$second_category<-substr(b[j],14,27)       
    q90_var_df_T2m$first_category<-first_category_name[i] 
    merge_q90_T2m<-rbind(merge_q90_T2m,q90_var_df_T2m)
    
    names(q90_var_df_DPT)<-c(1:60)
    q90_var_df_DPT$second_category<-substr(b[j],14,27)       
    q90_var_df_DPT$first_category<-first_category_name[i] 
    merge_q90_DPT<-rbind(merge_q90_DPT,q90_var_df_DPT)
    
    #new row for min
    names(min_var_df_WSPD)<-c(1:60)
    min_var_df_WSPD$second_category<-substr(b[j],14,27)       
    min_var_df_WSPD$first_category<-first_category_name[i] 
    merge_min_WSPD<-rbind(merge_min_WSPD,min_var_df_WSPD)
    
    names(min_var_df_T2m)<-c(1:60)
    min_var_df_T2m$second_category<-substr(b[j],14,27)       
    min_var_df_T2m$first_category<-first_category_name[i] 
    merge_min_T2m<-rbind(merge_min_T2m,min_var_df_T2m)
    
    names(min_var_df_DPT)<-c(1:60)
    min_var_df_DPT$second_category<-substr(b[j],14,27)       
    min_var_df_DPT$first_category<-first_category_name[i] 
    merge_min_DPT<-rbind(merge_min_DPT,min_var_df_DPT)
    
    #new row for max
    names(max_var_df_WSPD)<-c(1:60)
    max_var_df_WSPD$second_category<-substr(b[j],14,27)       
    max_var_df_WSPD$first_category<-first_category_name[i] 
    merge_max_WSPD<-rbind(merge_max_WSPD,max_var_df_WSPD)
    
    names(max_var_df_T2m)<-c(1:60)
    max_var_df_T2m$second_category<-substr(b[j],14,27)       
    max_var_df_T2m$first_category<-first_category_name[i] 
    merge_max_T2m<-rbind(merge_max_T2m,max_var_df_T2m)
    
    names(max_var_df_DPT)<-c(1:60)
    max_var_df_DPT$second_category<-substr(b[j],14,27)       
    max_var_df_DPT$first_category<-first_category_name[i] 
    merge_max_DPT<-rbind(merge_max_DPT,max_var_df_DPT)
    
  }
}   
Sys.time()
setwd("E:/ECMWFdata/pre-process")
write.table(merge_HRES_WSPD, "WSPD_6h_HRES_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_avg_WSPD, "WSPD_6h_mean_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_sd_WSPD, "WSPD_6h_sd_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_median_WSPD, "WSPD_6h_median_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q10_WSPD, "WSPD_6h_q10_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q25_WSPD, "WSPD_6h_q25_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q75_WSPD, "WSPD_6h_q75_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q90_WSPD, "WSPD_6h_q90_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_min_WSPD, "WSPD_6h_min_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_max_WSPD, "WSPD_6h_max_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)

write.table(merge_HRES_T2m, "T2m_6h_HRES_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_avg_T2m, "T2m_6h_mean_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_sd_T2m, "T2m_6h_sd_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_median_T2m, "T2m_6h_median_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q10_T2m, "T2m_6h_q10_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q25_T2m, "T2m_6h_q25_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q75_T2m, "T2m_6h_q75_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q90_T2m, "T2m_6h_q90_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_min_T2m, "T2m_6h_min_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_max_T2m, "T2m_6h_max_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)

write.table(merge_HRES_DPT, "DPT_6h_HRES_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_avg_DPT, "DPT_6h_mean_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_sd_DPT, "DPT_6h_sd_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_median_DPT, "DPT_6h_median_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q10_DPT, "DPT_6h_q10_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q25_DPT, "DPT_6h_q25_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q75_DPT, "DPT_6h_q75_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q90_DPT, "DPT_6h_q90_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_min_DPT, "DPT_6h_min_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_max_DPT, "DPT_6h_max_2021_Apr.txt", sep = "\t", col.names = T, row.names = F)

rm(list = ls())



#WSPD (m/s)
Sys.time()
setwd('E:/ECMWFdata/ECME_2021')
first_category_name = list.files("Jan")
dir = paste("./Jan/",first_category_name,sep="")
n = length(dir) 
n_sub<-rep(0,n)
n_sub<-as.data.frame(n_sub)
n_sub<-t(n_sub)
head(n_sub)
#create the df for HRES
merge_HRES_WSPD = as.data.frame(matrix(1:40,nrow = 1, ncol = 40))
names(merge_HRES_WSPD) = c(1:40)
merge_HRES_WSPD$second_category <- 'second_category'
merge_HRES_WSPD$first_category <- 'first_category'

merge_HRES_T2m = as.data.frame(matrix(1:40,nrow = 1, ncol = 40))
names(merge_HRES_T2m) = c(1:40)
merge_HRES_T2m$second_category <- 'second_category'
merge_HRES_T2m$first_category <- 'first_category'

merge_HRES_DPT = as.data.frame(matrix(1:40,nrow = 1, ncol = 40))
names(merge_HRES_DPT) = c(1:40)
merge_HRES_DPT$second_category <- 'second_category'
merge_HRES_DPT$first_category <- 'first_category'

#create the df for avg
merge_avg_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_avg_WSPD) = c(1:60)
merge_avg_WSPD$second_category <- 'second_category'
merge_avg_WSPD$first_category <- 'first_category'

merge_avg_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_avg_T2m) = c(1:60)
merge_avg_T2m$second_category <- 'second_category'
merge_avg_T2m$first_category <- 'first_category'

merge_avg_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_avg_DPT) = c(1:60)
merge_avg_DPT$second_category <- 'second_category'
merge_avg_DPT$first_category <- 'first_category'


#create the df for sd
merge_sd_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_sd_WSPD) = c(1:60)
merge_sd_WSPD$second_category <- 'second_category'
merge_sd_WSPD$first_category <- 'first_category'

merge_sd_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_sd_T2m) = c(1:60)
merge_sd_T2m$second_category <- 'second_category'
merge_sd_T2m$first_category <- 'first_category'

merge_sd_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_sd_DPT) = c(1:60)
merge_sd_DPT$second_category <- 'second_category'
merge_sd_DPT$first_category <- 'first_category'

#create the df for median
merge_median_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_median_WSPD) = c(1:60)
merge_median_WSPD$second_category <- 'second_category'
merge_median_WSPD$first_category <- 'first_category'

merge_median_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_median_T2m) = c(1:60)
merge_median_T2m$second_category <- 'second_category'
merge_median_T2m$first_category <- 'first_category'

merge_median_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_median_DPT) = c(1:60)
merge_median_DPT$second_category <- 'second_category'
merge_median_DPT$first_category <- 'first_category'

#create the df for q10
merge_q10_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q10_WSPD) = c(1:60)
merge_q10_WSPD$second_category <- 'second_category'
merge_q10_WSPD$first_category <- 'first_category'

merge_q10_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q10_T2m) = c(1:60)
merge_q10_T2m$second_category <- 'second_category'
merge_q10_T2m$first_category <- 'first_category'

merge_q10_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q10_DPT) = c(1:60)
merge_q10_DPT$second_category <- 'second_category'
merge_q10_DPT$first_category <- 'first_category'

#create the df for q25
merge_q25_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q25_WSPD) = c(1:60)
merge_q25_WSPD$second_category <- 'second_category'
merge_q25_WSPD$first_category <- 'first_category'

merge_q25_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q25_T2m) = c(1:60)
merge_q25_T2m$second_category <- 'second_category'
merge_q25_T2m$first_category <- 'first_category'

merge_q25_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q25_DPT) = c(1:60)
merge_q25_DPT$second_category <- 'second_category'
merge_q25_DPT$first_category <- 'first_category'

#createe the df for q75
merge_q75_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q75_WSPD) = c(1:60)
merge_q75_WSPD$second_category <- 'second_category'
merge_q75_WSPD$first_category <- 'first_category'

merge_q75_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q75_T2m) = c(1:60)
merge_q75_T2m$second_category <- 'second_category'
merge_q75_T2m$first_category <- 'first_category'

merge_q75_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q75_DPT) = c(1:60)
merge_q75_DPT$second_category <- 'second_category'
merge_q75_DPT$first_category <- 'first_category'


#create the df for q90
merge_q90_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q90_WSPD) = c(1:60)
merge_q90_WSPD$second_category <- 'second_category'
merge_q90_WSPD$first_category <- 'first_category'

merge_q90_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q90_T2m) = c(1:60)
merge_q90_T2m$second_category <- 'second_category'
merge_q90_T2m$first_category <- 'first_category'

merge_q90_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q90_DPT) = c(1:60)
merge_q90_DPT$second_category <- 'second_category'
merge_q90_DPT$first_category <- 'first_category'


#create the df for min
merge_min_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_min_WSPD) = c(1:60)
merge_min_WSPD$second_category <- 'second_category'
merge_min_WSPD$first_category <- 'first_category'

merge_min_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_min_T2m) = c(1:60)
merge_min_T2m$second_category <- 'second_category'
merge_min_T2m$first_category <- 'first_category'

merge_min_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_min_DPT) = c(1:60)
merge_min_DPT$second_category <- 'second_category'
merge_min_DPT$first_category <- 'first_category'


#create the df for max
merge_max_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_max_WSPD) = c(1:60)
merge_max_WSPD$second_category <- 'second_category'
merge_max_WSPD$first_category <- 'first_category'

merge_max_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_max_T2m) = c(1:60)
merge_max_T2m$second_category <- 'second_category'
merge_max_T2m$first_category <- 'first_category'

merge_max_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_max_DPT) = c(1:60)
merge_max_DPT$second_category <- 'second_category'
merge_max_DPT$first_category <- 'first_category'






for(i in 1:n){         #对于每个一级目录(文件夹)
  b=list.files(dir[i]) #b是列出每个一级目录(文件夹)中每个xlsx文件的名称
  n_sub[i]=length(b)   #得到一级目录(文件夹)下xlsx的文件个数:n_sub
  
  for(j in 1:n_sub[i]){     #对于每个一级目录(文件夹)下的每个xlsx文件
    data<-read.delim(file=paste(dir[i],'/',b[j],sep=''),sep = "", header = F) #读取xlsx文件
    data_WSPD = data[161:213,]
    data_WSPD[,48:83]=trunc(data_WSPD[,48:83]/1000000)
    data_WSPD[2,64]=data_WSPD[2,64]/10^234
    data_WSPD = as.data.frame(lapply(data_WSPD,as.numeric))
    data_WSPD[2:53,]=data_WSPD[2:53,]/10
    WSPD_6hstep = data_WSPD[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
    #T2m (degrees Celcius)
    data_T2m = data[267:319,]
    data_T2m[,48:83]=trunc(data_T2m[,48:83]/1000000)
    data_T2m[2,64]=data_T2m[2,64]/10^234
    data_T2m = as.data.frame(lapply(data_T2m,as.numeric))
    data_T2m[2:53,]=data_T2m[2:53,]/10
    T2m_6hstep = data_T2m[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
    #DPT (degrees Celcius)
    data_DPT = data[320:372,]
    data_DPT[,48:83]=trunc(data_DPT[,48:83]/1000000)
    data_DPT[2,64]=data_DPT[2,64]/10^234
    data_DPT = as.data.frame(lapply(data_DPT,as.numeric))
    data_DPT[2:53,]=data_DPT[2:53,]/10
    DPT_6hstep = data_DPT[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
    
    WSPD_HRES = WSPD_6hstep[2,1:40]
    DPT_HRES = DPT_6hstep[2,1:40]
    T2m_HRES = T2m_6hstep[2,1:40]
    
    WSPD_ENS = WSPD_6hstep[3:53,]
    WSPD_ENS = as.data.frame(lapply(WSPD_ENS,as.numeric))
    
    T2m_ENS = T2m_6hstep[3:53,]
    T2m_ENS = as.data.frame(lapply(T2m_ENS,as.numeric))
    
    DPT_ENS = DPT_6hstep[3:53,]
    DPT_ENS = as.data.frame(lapply(DPT_ENS,as.numeric))
    
    avg_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    avg_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    avg_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    sd_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    sd_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    sd_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    
    median_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    median_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    median_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    
    q10_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q10_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q10_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    q25_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q25_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q25_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    q75_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q75_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q75_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    
    q90_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q90_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q90_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    min_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    min_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    min_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    max_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    max_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    max_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    
    
    for(k in 1:60){ 
      avg_var_df_WSPD[,k] <- mean(WSPD_ENS[, k])
      avg_var_df_T2m[,k] <- mean(T2m_ENS[, k])
      avg_var_df_DPT[,k] <- mean(DPT_ENS[, k])
      
      sd_var_df_WSPD[,k] <- sd(WSPD_ENS[, k])
      sd_var_df_T2m[,k] <- sd(T2m_ENS[, k])
      sd_var_df_DPT[,k] <- sd(DPT_ENS[, k])
      
      median_var_df_WSPD[,k] <- median(WSPD_ENS[, k])
      median_var_df_T2m[,k] <- median(T2m_ENS[, k])
      median_var_df_DPT[,k] <- median(DPT_ENS[, k])
      
      q10_var_df_WSPD[,k] <- quantile(WSPD_ENS[,k],probs=0.1)
      q10_var_df_T2m[,k] <- quantile(T2m_ENS[,k],probs=0.1)
      q10_var_df_DPT[,k] <- quantile(DPT_ENS[,k],probs=0.1)
      
      q25_var_df_WSPD[,k] <- quantile(WSPD_ENS[,k],probs=0.25)
      q25_var_df_T2m[,k] <- quantile(T2m_ENS[,k],probs=0.25)
      q25_var_df_DPT[,k] <- quantile(DPT_ENS[,k],probs=0.25)
      
      q75_var_df_WSPD[,k] <- quantile(WSPD_ENS[,k],probs=0.75)
      q75_var_df_T2m[,k] <- quantile(T2m_ENS[,k],probs=0.75)
      q75_var_df_DPT[,k] <- quantile(DPT_ENS[,k],probs=0.75)
      
      q90_var_df_WSPD[,k] <- quantile(WSPD_ENS[,k],probs=0.9)
      q90_var_df_T2m[,k] <- quantile(T2m_ENS[,k],probs=0.9)
      q90_var_df_DPT[,k] <- quantile(DPT_ENS[,k],probs=0.9)
      
      
      min_var_df_WSPD[,k] <- min(WSPD_ENS[, k])
      min_var_df_T2m[,k] <- min(T2m_ENS[, k])
      min_var_df_DPT[,k] <- min(DPT_ENS[, k])
      
      max_var_df_WSPD[,k] <- max(WSPD_ENS[, k])
      max_var_df_T2m[,k] <- max(T2m_ENS[, k])
      max_var_df_DPT[,k] <- max(DPT_ENS[, k])
      
    }    
    #new row for HRES
    names(WSPD_HRES)<-c(1:40)
    WSPD_HRES$second_category<-substr(b[j],14,27)      
    WSPD_HRES$first_category<-first_category_name[i]   
    merge_HRES_WSPD<-rbind(merge_HRES_WSPD,WSPD_HRES)
    
    names(T2m_HRES)<-c(1:40)
    T2m_HRES$second_category<-substr(b[j],14,27)      
    T2m_HRES$first_category<-first_category_name[i]   
    merge_HRES_T2m<-rbind(merge_HRES_T2m,T2m_HRES)
    
    names(DPT_HRES)<-c(1:40)
    DPT_HRES$second_category<-substr(b[j],14,27)      
    DPT_HRES$first_category<-first_category_name[i]   
    merge_HRES_DPT<-rbind(merge_HRES_DPT,DPT_HRES)
    
    #new row for avg
    names(avg_var_df_WSPD)<-c(1:60)
    avg_var_df_WSPD$second_category<-substr(b[j],14,27)       
    avg_var_df_WSPD$first_category<-first_category_name[i] 
    merge_avg_WSPD<-rbind(merge_avg_WSPD,avg_var_df_WSPD)
    
    names(avg_var_df_T2m)<-c(1:60)
    avg_var_df_T2m$second_category<-substr(b[j],14,27)       
    avg_var_df_T2m$first_category<-first_category_name[i] 
    merge_avg_T2m<-rbind(merge_avg_T2m,avg_var_df_T2m)
    
    names(avg_var_df_DPT)<-c(1:60)
    avg_var_df_DPT$second_category<-substr(b[j],14,27)       
    avg_var_df_DPT$first_category<-first_category_name[i] 
    merge_avg_DPT<-rbind(merge_avg_DPT,avg_var_df_DPT)
    
    #new row for sd
    names(sd_var_df_WSPD)<-c(1:60)
    sd_var_df_WSPD$second_category<-substr(b[j],14,27)       
    sd_var_df_WSPD$first_category<-first_category_name[i] 
    merge_sd_WSPD<-rbind(merge_sd_WSPD,sd_var_df_WSPD)
    
    names(sd_var_df_T2m)<-c(1:60)
    sd_var_df_T2m$second_category<-substr(b[j],14,27)       
    sd_var_df_T2m$first_category<-first_category_name[i] 
    merge_sd_T2m<-rbind(merge_sd_T2m,sd_var_df_T2m)
    
    names(sd_var_df_DPT)<-c(1:60)
    sd_var_df_DPT$second_category<-substr(b[j],14,27)       
    sd_var_df_DPT$first_category<-first_category_name[i] 
    merge_sd_DPT<-rbind(merge_sd_DPT,sd_var_df_DPT)
    
    
    #new row for median
    names(median_var_df_WSPD)<-c(1:60)
    median_var_df_WSPD$second_category<-substr(b[j],14,27)       
    median_var_df_WSPD$first_category<-first_category_name[i] 
    merge_median_WSPD<-rbind(merge_median_WSPD,median_var_df_WSPD)
    
    names(median_var_df_T2m)<-c(1:60)
    median_var_df_T2m$second_category<-substr(b[j],14,27)       
    median_var_df_T2m$first_category<-first_category_name[i] 
    merge_median_T2m<-rbind(merge_median_T2m,median_var_df_T2m)
    
    names(median_var_df_DPT)<-c(1:60)
    median_var_df_DPT$second_category<-substr(b[j],14,27)       
    median_var_df_DPT$first_category<-first_category_name[i] 
    merge_median_DPT<-rbind(merge_median_DPT,median_var_df_DPT)
    
    #new row for q10
    names(q10_var_df_WSPD)<-c(1:60)
    q10_var_df_WSPD$second_category<-substr(b[j],14,27)       
    q10_var_df_WSPD$first_category<-first_category_name[i] 
    merge_q10_WSPD<-rbind(merge_q10_WSPD,q10_var_df_WSPD)
    
    names(q10_var_df_T2m)<-c(1:60)
    q10_var_df_T2m$second_category<-substr(b[j],14,27)       
    q10_var_df_T2m$first_category<-first_category_name[i] 
    merge_q10_T2m<-rbind(merge_q10_T2m,q10_var_df_T2m)
    
    names(q10_var_df_DPT)<-c(1:60)
    q10_var_df_DPT$second_category<-substr(b[j],14,27)       
    q10_var_df_DPT$first_category<-first_category_name[i] 
    merge_q10_DPT<-rbind(merge_q10_DPT,q10_var_df_DPT)
    
    #new row for q25
    names(q25_var_df_WSPD)<-c(1:60)
    q25_var_df_WSPD$second_category<-substr(b[j],14,27)       
    q25_var_df_WSPD$first_category<-first_category_name[i] 
    merge_q25_WSPD<-rbind(merge_q25_WSPD,q25_var_df_WSPD)
    
    names(q25_var_df_T2m)<-c(1:60)
    q25_var_df_T2m$second_category<-substr(b[j],14,27)       
    q25_var_df_T2m$first_category<-first_category_name[i] 
    merge_q25_T2m<-rbind(merge_q25_T2m,q25_var_df_T2m)
    
    names(q25_var_df_DPT)<-c(1:60)
    q25_var_df_DPT$second_category<-substr(b[j],14,27)       
    q25_var_df_DPT$first_category<-first_category_name[i] 
    merge_q25_DPT<-rbind(merge_q25_DPT,q25_var_df_DPT)
    
    #new row for q75
    names(q75_var_df_WSPD)<-c(1:60)
    q75_var_df_WSPD$second_category<-substr(b[j],14,27)       
    q75_var_df_WSPD$first_category<-first_category_name[i] 
    merge_q75_WSPD<-rbind(merge_q75_WSPD,q75_var_df_WSPD)
    
    names(q75_var_df_T2m)<-c(1:60)
    q75_var_df_T2m$second_category<-substr(b[j],14,27)       
    q75_var_df_T2m$first_category<-first_category_name[i] 
    merge_q75_T2m<-rbind(merge_q75_T2m,q75_var_df_T2m)
    
    names(q75_var_df_DPT)<-c(1:60)
    q75_var_df_DPT$second_category<-substr(b[j],14,27)       
    q75_var_df_DPT$first_category<-first_category_name[i] 
    merge_q75_DPT<-rbind(merge_q75_DPT,q75_var_df_DPT)
    
    
    
    #new row for q90
    names(q90_var_df_WSPD)<-c(1:60)
    q90_var_df_WSPD$second_category<-substr(b[j],14,27)       
    q90_var_df_WSPD$first_category<-first_category_name[i] 
    merge_q90_WSPD<-rbind(merge_q90_WSPD,q90_var_df_WSPD)
    
    names(q90_var_df_T2m)<-c(1:60)
    q90_var_df_T2m$second_category<-substr(b[j],14,27)       
    q90_var_df_T2m$first_category<-first_category_name[i] 
    merge_q90_T2m<-rbind(merge_q90_T2m,q90_var_df_T2m)
    
    names(q90_var_df_DPT)<-c(1:60)
    q90_var_df_DPT$second_category<-substr(b[j],14,27)       
    q90_var_df_DPT$first_category<-first_category_name[i] 
    merge_q90_DPT<-rbind(merge_q90_DPT,q90_var_df_DPT)
    
    #new row for min
    names(min_var_df_WSPD)<-c(1:60)
    min_var_df_WSPD$second_category<-substr(b[j],14,27)       
    min_var_df_WSPD$first_category<-first_category_name[i] 
    merge_min_WSPD<-rbind(merge_min_WSPD,min_var_df_WSPD)
    
    names(min_var_df_T2m)<-c(1:60)
    min_var_df_T2m$second_category<-substr(b[j],14,27)       
    min_var_df_T2m$first_category<-first_category_name[i] 
    merge_min_T2m<-rbind(merge_min_T2m,min_var_df_T2m)
    
    names(min_var_df_DPT)<-c(1:60)
    min_var_df_DPT$second_category<-substr(b[j],14,27)       
    min_var_df_DPT$first_category<-first_category_name[i] 
    merge_min_DPT<-rbind(merge_min_DPT,min_var_df_DPT)
    
    #new row for max
    names(max_var_df_WSPD)<-c(1:60)
    max_var_df_WSPD$second_category<-substr(b[j],14,27)       
    max_var_df_WSPD$first_category<-first_category_name[i] 
    merge_max_WSPD<-rbind(merge_max_WSPD,max_var_df_WSPD)
    
    names(max_var_df_T2m)<-c(1:60)
    max_var_df_T2m$second_category<-substr(b[j],14,27)       
    max_var_df_T2m$first_category<-first_category_name[i] 
    merge_max_T2m<-rbind(merge_max_T2m,max_var_df_T2m)
    
    names(max_var_df_DPT)<-c(1:60)
    max_var_df_DPT$second_category<-substr(b[j],14,27)       
    max_var_df_DPT$first_category<-first_category_name[i] 
    merge_max_DPT<-rbind(merge_max_DPT,max_var_df_DPT)
    
  }
}   
Sys.time()
setwd("E:/ECMWFdata/pre-process")
write.table(merge_HRES_WSPD, "WSPD_6h_HRES_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_avg_WSPD, "WSPD_6h_mean_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_sd_WSPD, "WSPD_6h_sd_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_median_WSPD, "WSPD_6h_median_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q10_WSPD, "WSPD_6h_q10_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q25_WSPD, "WSPD_6h_q25_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q75_WSPD, "WSPD_6h_q75_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q90_WSPD, "WSPD_6h_q90_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_min_WSPD, "WSPD_6h_min_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_max_WSPD, "WSPD_6h_max_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)

write.table(merge_HRES_T2m, "T2m_6h_HRES_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_avg_T2m, "T2m_6h_mean_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_sd_T2m, "T2m_6h_sd_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_median_T2m, "T2m_6h_median_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q10_T2m, "T2m_6h_q10_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q25_T2m, "T2m_6h_q25_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q75_T2m, "T2m_6h_q75_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q90_T2m, "T2m_6h_q90_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_min_T2m, "T2m_6h_min_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_max_T2m, "T2m_6h_max_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)

write.table(merge_HRES_DPT, "DPT_6h_HRES_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_avg_DPT, "DPT_6h_mean_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_sd_DPT, "DPT_6h_sd_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_median_DPT, "DPT_6h_median_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q10_DPT, "DPT_6h_q10_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q25_DPT, "DPT_6h_q25_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q75_DPT, "DPT_6h_q75_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q90_DPT, "DPT_6h_q90_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_min_DPT, "DPT_6h_min_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_max_DPT, "DPT_6h_max_2021_Jan.txt", sep = "\t", col.names = T, row.names = F)


rm(list = ls())



#WSPD (m/s)
Sys.time()
setwd('E:/ECMWFdata/ECME_2021')
first_category_name = list.files("Feb")
dir = paste("./Feb/",first_category_name,sep="")
n = length(dir) 
n_sub<-rep(0,n)
n_sub<-as.data.frame(n_sub)
n_sub<-t(n_sub)
head(n_sub)
#create the df for HRES
merge_HRES_WSPD = as.data.frame(matrix(1:40,nrow = 1, ncol = 40))
names(merge_HRES_WSPD) = c(1:40)
merge_HRES_WSPD$second_category <- 'second_category'
merge_HRES_WSPD$first_category <- 'first_category'

merge_HRES_T2m = as.data.frame(matrix(1:40,nrow = 1, ncol = 40))
names(merge_HRES_T2m) = c(1:40)
merge_HRES_T2m$second_category <- 'second_category'
merge_HRES_T2m$first_category <- 'first_category'

merge_HRES_DPT = as.data.frame(matrix(1:40,nrow = 1, ncol = 40))
names(merge_HRES_DPT) = c(1:40)
merge_HRES_DPT$second_category <- 'second_category'
merge_HRES_DPT$first_category <- 'first_category'

#create the df for avg
merge_avg_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_avg_WSPD) = c(1:60)
merge_avg_WSPD$second_category <- 'second_category'
merge_avg_WSPD$first_category <- 'first_category'

merge_avg_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_avg_T2m) = c(1:60)
merge_avg_T2m$second_category <- 'second_category'
merge_avg_T2m$first_category <- 'first_category'

merge_avg_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_avg_DPT) = c(1:60)
merge_avg_DPT$second_category <- 'second_category'
merge_avg_DPT$first_category <- 'first_category'


#create the df for sd
merge_sd_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_sd_WSPD) = c(1:60)
merge_sd_WSPD$second_category <- 'second_category'
merge_sd_WSPD$first_category <- 'first_category'

merge_sd_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_sd_T2m) = c(1:60)
merge_sd_T2m$second_category <- 'second_category'
merge_sd_T2m$first_category <- 'first_category'

merge_sd_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_sd_DPT) = c(1:60)
merge_sd_DPT$second_category <- 'second_category'
merge_sd_DPT$first_category <- 'first_category'

#create the df for median
merge_median_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_median_WSPD) = c(1:60)
merge_median_WSPD$second_category <- 'second_category'
merge_median_WSPD$first_category <- 'first_category'

merge_median_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_median_T2m) = c(1:60)
merge_median_T2m$second_category <- 'second_category'
merge_median_T2m$first_category <- 'first_category'

merge_median_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_median_DPT) = c(1:60)
merge_median_DPT$second_category <- 'second_category'
merge_median_DPT$first_category <- 'first_category'

#create the df for q10
merge_q10_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q10_WSPD) = c(1:60)
merge_q10_WSPD$second_category <- 'second_category'
merge_q10_WSPD$first_category <- 'first_category'

merge_q10_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q10_T2m) = c(1:60)
merge_q10_T2m$second_category <- 'second_category'
merge_q10_T2m$first_category <- 'first_category'

merge_q10_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q10_DPT) = c(1:60)
merge_q10_DPT$second_category <- 'second_category'
merge_q10_DPT$first_category <- 'first_category'

#create the df for q25
merge_q25_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q25_WSPD) = c(1:60)
merge_q25_WSPD$second_category <- 'second_category'
merge_q25_WSPD$first_category <- 'first_category'

merge_q25_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q25_T2m) = c(1:60)
merge_q25_T2m$second_category <- 'second_category'
merge_q25_T2m$first_category <- 'first_category'

merge_q25_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q25_DPT) = c(1:60)
merge_q25_DPT$second_category <- 'second_category'
merge_q25_DPT$first_category <- 'first_category'

#createe the df for q75
merge_q75_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q75_WSPD) = c(1:60)
merge_q75_WSPD$second_category <- 'second_category'
merge_q75_WSPD$first_category <- 'first_category'

merge_q75_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q75_T2m) = c(1:60)
merge_q75_T2m$second_category <- 'second_category'
merge_q75_T2m$first_category <- 'first_category'

merge_q75_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q75_DPT) = c(1:60)
merge_q75_DPT$second_category <- 'second_category'
merge_q75_DPT$first_category <- 'first_category'


#create the df for q90
merge_q90_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q90_WSPD) = c(1:60)
merge_q90_WSPD$second_category <- 'second_category'
merge_q90_WSPD$first_category <- 'first_category'

merge_q90_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q90_T2m) = c(1:60)
merge_q90_T2m$second_category <- 'second_category'
merge_q90_T2m$first_category <- 'first_category'

merge_q90_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q90_DPT) = c(1:60)
merge_q90_DPT$second_category <- 'second_category'
merge_q90_DPT$first_category <- 'first_category'


#create the df for min
merge_min_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_min_WSPD) = c(1:60)
merge_min_WSPD$second_category <- 'second_category'
merge_min_WSPD$first_category <- 'first_category'

merge_min_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_min_T2m) = c(1:60)
merge_min_T2m$second_category <- 'second_category'
merge_min_T2m$first_category <- 'first_category'

merge_min_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_min_DPT) = c(1:60)
merge_min_DPT$second_category <- 'second_category'
merge_min_DPT$first_category <- 'first_category'


#create the df for max
merge_max_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_max_WSPD) = c(1:60)
merge_max_WSPD$second_category <- 'second_category'
merge_max_WSPD$first_category <- 'first_category'

merge_max_T2m = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_max_T2m) = c(1:60)
merge_max_T2m$second_category <- 'second_category'
merge_max_T2m$first_category <- 'first_category'

merge_max_DPT = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_max_DPT) = c(1:60)
merge_max_DPT$second_category <- 'second_category'
merge_max_DPT$first_category <- 'first_category'






for(i in 1:n){         #对于每个一级目录(文件夹)
  b=list.files(dir[i]) #b是列出每个一级目录(文件夹)中每个xlsx文件的名称
  n_sub[i]=length(b)   #得到一级目录(文件夹)下xlsx的文件个数:n_sub
  
  for(j in 1:n_sub[i]){     #对于每个一级目录(文件夹)下的每个xlsx文件
    data<-read.delim(file=paste(dir[i],'/',b[j],sep=''),sep = "", header = F) #读取xlsx文件
    data_WSPD = data[161:213,]
    data_WSPD[,48:83]=trunc(data_WSPD[,48:83]/1000000)
    data_WSPD[2,64]=data_WSPD[2,64]/10^234
    data_WSPD = as.data.frame(lapply(data_WSPD,as.numeric))
    data_WSPD[2:53,]=data_WSPD[2:53,]/10
    WSPD_6hstep = data_WSPD[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
    #T2m (degrees Celcius)
    data_T2m = data[267:319,]
    data_T2m[,48:83]=trunc(data_T2m[,48:83]/1000000)
    data_T2m[2,64]=data_T2m[2,64]/10^234
    data_T2m = as.data.frame(lapply(data_T2m,as.numeric))
    data_T2m[2:53,]=data_T2m[2:53,]/10
    T2m_6hstep = data_T2m[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
    #DPT (degrees Celcius)
    data_DPT = data[320:372,]
    data_DPT[,48:83]=trunc(data_DPT[,48:83]/1000000)
    data_DPT[2,64]=data_DPT[2,64]/10^234
    data_DPT = as.data.frame(lapply(data_DPT,as.numeric))
    data_DPT[2:53,]=data_DPT[2:53,]/10
    DPT_6hstep = data_DPT[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
    
    WSPD_HRES = WSPD_6hstep[2,1:40]
    DPT_HRES = DPT_6hstep[2,1:40]
    T2m_HRES = T2m_6hstep[2,1:40]
    
    WSPD_ENS = WSPD_6hstep[3:53,]
    WSPD_ENS = as.data.frame(lapply(WSPD_ENS,as.numeric))
    
    T2m_ENS = T2m_6hstep[3:53,]
    T2m_ENS = as.data.frame(lapply(T2m_ENS,as.numeric))
    
    DPT_ENS = DPT_6hstep[3:53,]
    DPT_ENS = as.data.frame(lapply(DPT_ENS,as.numeric))
    
    avg_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    avg_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    avg_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    sd_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    sd_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    sd_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    
    median_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    median_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    median_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    
    q10_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q10_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q10_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    q25_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q25_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q25_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    q75_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q75_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q75_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    
    q90_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q90_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    q90_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    min_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    min_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    min_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    max_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    max_var_df_T2m <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    max_var_df_DPT <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
    
    
    
    for(k in 1:60){ 
      avg_var_df_WSPD[,k] <- mean(WSPD_ENS[, k])
      avg_var_df_T2m[,k] <- mean(T2m_ENS[, k])
      avg_var_df_DPT[,k] <- mean(DPT_ENS[, k])
      
      sd_var_df_WSPD[,k] <- sd(WSPD_ENS[, k])
      sd_var_df_T2m[,k] <- sd(T2m_ENS[, k])
      sd_var_df_DPT[,k] <- sd(DPT_ENS[, k])
      
      median_var_df_WSPD[,k] <- median(WSPD_ENS[, k])
      median_var_df_T2m[,k] <- median(T2m_ENS[, k])
      median_var_df_DPT[,k] <- median(DPT_ENS[, k])
      
      q10_var_df_WSPD[,k] <- quantile(WSPD_ENS[,k],probs=0.1)
      q10_var_df_T2m[,k] <- quantile(T2m_ENS[,k],probs=0.1)
      q10_var_df_DPT[,k] <- quantile(DPT_ENS[,k],probs=0.1)
      
      q25_var_df_WSPD[,k] <- quantile(WSPD_ENS[,k],probs=0.25)
      q25_var_df_T2m[,k] <- quantile(T2m_ENS[,k],probs=0.25)
      q25_var_df_DPT[,k] <- quantile(DPT_ENS[,k],probs=0.25)
      
      q75_var_df_WSPD[,k] <- quantile(WSPD_ENS[,k],probs=0.75)
      q75_var_df_T2m[,k] <- quantile(T2m_ENS[,k],probs=0.75)
      q75_var_df_DPT[,k] <- quantile(DPT_ENS[,k],probs=0.75)
      
      q90_var_df_WSPD[,k] <- quantile(WSPD_ENS[,k],probs=0.9)
      q90_var_df_T2m[,k] <- quantile(T2m_ENS[,k],probs=0.9)
      q90_var_df_DPT[,k] <- quantile(DPT_ENS[,k],probs=0.9)
      
      
      min_var_df_WSPD[,k] <- min(WSPD_ENS[, k])
      min_var_df_T2m[,k] <- min(T2m_ENS[, k])
      min_var_df_DPT[,k] <- min(DPT_ENS[, k])
      
      max_var_df_WSPD[,k] <- max(WSPD_ENS[, k])
      max_var_df_T2m[,k] <- max(T2m_ENS[, k])
      max_var_df_DPT[,k] <- max(DPT_ENS[, k])
      
    }    
    #new row for HRES
    names(WSPD_HRES)<-c(1:40)
    WSPD_HRES$second_category<-substr(b[j],14,27)      
    WSPD_HRES$first_category<-first_category_name[i]   
    merge_HRES_WSPD<-rbind(merge_HRES_WSPD,WSPD_HRES)
    
    names(T2m_HRES)<-c(1:40)
    T2m_HRES$second_category<-substr(b[j],14,27)      
    T2m_HRES$first_category<-first_category_name[i]   
    merge_HRES_T2m<-rbind(merge_HRES_T2m,T2m_HRES)
    
    names(DPT_HRES)<-c(1:40)
    DPT_HRES$second_category<-substr(b[j],14,27)      
    DPT_HRES$first_category<-first_category_name[i]   
    merge_HRES_DPT<-rbind(merge_HRES_DPT,DPT_HRES)
    
    #new row for avg
    names(avg_var_df_WSPD)<-c(1:60)
    avg_var_df_WSPD$second_category<-substr(b[j],14,27)       
    avg_var_df_WSPD$first_category<-first_category_name[i] 
    merge_avg_WSPD<-rbind(merge_avg_WSPD,avg_var_df_WSPD)
    
    names(avg_var_df_T2m)<-c(1:60)
    avg_var_df_T2m$second_category<-substr(b[j],14,27)       
    avg_var_df_T2m$first_category<-first_category_name[i] 
    merge_avg_T2m<-rbind(merge_avg_T2m,avg_var_df_T2m)
    
    names(avg_var_df_DPT)<-c(1:60)
    avg_var_df_DPT$second_category<-substr(b[j],14,27)       
    avg_var_df_DPT$first_category<-first_category_name[i] 
    merge_avg_DPT<-rbind(merge_avg_DPT,avg_var_df_DPT)
    
    #new row for sd
    names(sd_var_df_WSPD)<-c(1:60)
    sd_var_df_WSPD$second_category<-substr(b[j],14,27)       
    sd_var_df_WSPD$first_category<-first_category_name[i] 
    merge_sd_WSPD<-rbind(merge_sd_WSPD,sd_var_df_WSPD)
    
    names(sd_var_df_T2m)<-c(1:60)
    sd_var_df_T2m$second_category<-substr(b[j],14,27)       
    sd_var_df_T2m$first_category<-first_category_name[i] 
    merge_sd_T2m<-rbind(merge_sd_T2m,sd_var_df_T2m)
    
    names(sd_var_df_DPT)<-c(1:60)
    sd_var_df_DPT$second_category<-substr(b[j],14,27)       
    sd_var_df_DPT$first_category<-first_category_name[i] 
    merge_sd_DPT<-rbind(merge_sd_DPT,sd_var_df_DPT)
    
    
    #new row for median
    names(median_var_df_WSPD)<-c(1:60)
    median_var_df_WSPD$second_category<-substr(b[j],14,27)       
    median_var_df_WSPD$first_category<-first_category_name[i] 
    merge_median_WSPD<-rbind(merge_median_WSPD,median_var_df_WSPD)
    
    names(median_var_df_T2m)<-c(1:60)
    median_var_df_T2m$second_category<-substr(b[j],14,27)       
    median_var_df_T2m$first_category<-first_category_name[i] 
    merge_median_T2m<-rbind(merge_median_T2m,median_var_df_T2m)
    
    names(median_var_df_DPT)<-c(1:60)
    median_var_df_DPT$second_category<-substr(b[j],14,27)       
    median_var_df_DPT$first_category<-first_category_name[i] 
    merge_median_DPT<-rbind(merge_median_DPT,median_var_df_DPT)
    
    #new row for q10
    names(q10_var_df_WSPD)<-c(1:60)
    q10_var_df_WSPD$second_category<-substr(b[j],14,27)       
    q10_var_df_WSPD$first_category<-first_category_name[i] 
    merge_q10_WSPD<-rbind(merge_q10_WSPD,q10_var_df_WSPD)
    
    names(q10_var_df_T2m)<-c(1:60)
    q10_var_df_T2m$second_category<-substr(b[j],14,27)       
    q10_var_df_T2m$first_category<-first_category_name[i] 
    merge_q10_T2m<-rbind(merge_q10_T2m,q10_var_df_T2m)
    
    names(q10_var_df_DPT)<-c(1:60)
    q10_var_df_DPT$second_category<-substr(b[j],14,27)       
    q10_var_df_DPT$first_category<-first_category_name[i] 
    merge_q10_DPT<-rbind(merge_q10_DPT,q10_var_df_DPT)
    
    #new row for q25
    names(q25_var_df_WSPD)<-c(1:60)
    q25_var_df_WSPD$second_category<-substr(b[j],14,27)       
    q25_var_df_WSPD$first_category<-first_category_name[i] 
    merge_q25_WSPD<-rbind(merge_q25_WSPD,q25_var_df_WSPD)
    
    names(q25_var_df_T2m)<-c(1:60)
    q25_var_df_T2m$second_category<-substr(b[j],14,27)       
    q25_var_df_T2m$first_category<-first_category_name[i] 
    merge_q25_T2m<-rbind(merge_q25_T2m,q25_var_df_T2m)
    
    names(q25_var_df_DPT)<-c(1:60)
    q25_var_df_DPT$second_category<-substr(b[j],14,27)       
    q25_var_df_DPT$first_category<-first_category_name[i] 
    merge_q25_DPT<-rbind(merge_q25_DPT,q25_var_df_DPT)
    
    #new row for q75
    names(q75_var_df_WSPD)<-c(1:60)
    q75_var_df_WSPD$second_category<-substr(b[j],14,27)       
    q75_var_df_WSPD$first_category<-first_category_name[i] 
    merge_q75_WSPD<-rbind(merge_q75_WSPD,q75_var_df_WSPD)
    
    names(q75_var_df_T2m)<-c(1:60)
    q75_var_df_T2m$second_category<-substr(b[j],14,27)       
    q75_var_df_T2m$first_category<-first_category_name[i] 
    merge_q75_T2m<-rbind(merge_q75_T2m,q75_var_df_T2m)
    
    names(q75_var_df_DPT)<-c(1:60)
    q75_var_df_DPT$second_category<-substr(b[j],14,27)       
    q75_var_df_DPT$first_category<-first_category_name[i] 
    merge_q75_DPT<-rbind(merge_q75_DPT,q75_var_df_DPT)
    
    
    
    #new row for q90
    names(q90_var_df_WSPD)<-c(1:60)
    q90_var_df_WSPD$second_category<-substr(b[j],14,27)       
    q90_var_df_WSPD$first_category<-first_category_name[i] 
    merge_q90_WSPD<-rbind(merge_q90_WSPD,q90_var_df_WSPD)
    
    names(q90_var_df_T2m)<-c(1:60)
    q90_var_df_T2m$second_category<-substr(b[j],14,27)       
    q90_var_df_T2m$first_category<-first_category_name[i] 
    merge_q90_T2m<-rbind(merge_q90_T2m,q90_var_df_T2m)
    
    names(q90_var_df_DPT)<-c(1:60)
    q90_var_df_DPT$second_category<-substr(b[j],14,27)       
    q90_var_df_DPT$first_category<-first_category_name[i] 
    merge_q90_DPT<-rbind(merge_q90_DPT,q90_var_df_DPT)
    
    #new row for min
    names(min_var_df_WSPD)<-c(1:60)
    min_var_df_WSPD$second_category<-substr(b[j],14,27)       
    min_var_df_WSPD$first_category<-first_category_name[i] 
    merge_min_WSPD<-rbind(merge_min_WSPD,min_var_df_WSPD)
    
    names(min_var_df_T2m)<-c(1:60)
    min_var_df_T2m$second_category<-substr(b[j],14,27)       
    min_var_df_T2m$first_category<-first_category_name[i] 
    merge_min_T2m<-rbind(merge_min_T2m,min_var_df_T2m)
    
    names(min_var_df_DPT)<-c(1:60)
    min_var_df_DPT$second_category<-substr(b[j],14,27)       
    min_var_df_DPT$first_category<-first_category_name[i] 
    merge_min_DPT<-rbind(merge_min_DPT,min_var_df_DPT)
    
    #new row for max
    names(max_var_df_WSPD)<-c(1:60)
    max_var_df_WSPD$second_category<-substr(b[j],14,27)       
    max_var_df_WSPD$first_category<-first_category_name[i] 
    merge_max_WSPD<-rbind(merge_max_WSPD,max_var_df_WSPD)
    
    names(max_var_df_T2m)<-c(1:60)
    max_var_df_T2m$second_category<-substr(b[j],14,27)       
    max_var_df_T2m$first_category<-first_category_name[i] 
    merge_max_T2m<-rbind(merge_max_T2m,max_var_df_T2m)
    
    names(max_var_df_DPT)<-c(1:60)
    max_var_df_DPT$second_category<-substr(b[j],14,27)       
    max_var_df_DPT$first_category<-first_category_name[i] 
    merge_max_DPT<-rbind(merge_max_DPT,max_var_df_DPT)
    
  }
}   
Sys.time()
setwd("E:/ECMWFdata/pre-process")
write.table(merge_HRES_WSPD, "WSPD_6h_HRES_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_avg_WSPD, "WSPD_6h_mean_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_sd_WSPD, "WSPD_6h_sd_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_median_WSPD, "WSPD_6h_median_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q10_WSPD, "WSPD_6h_q10_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q25_WSPD, "WSPD_6h_q25_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q75_WSPD, "WSPD_6h_q75_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q90_WSPD, "WSPD_6h_q90_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_min_WSPD, "WSPD_6h_min_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_max_WSPD, "WSPD_6h_max_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)

write.table(merge_HRES_T2m, "T2m_6h_HRES_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_avg_T2m, "T2m_6h_mean_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_sd_T2m, "T2m_6h_sd_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_median_T2m, "T2m_6h_median_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q10_T2m, "T2m_6h_q10_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q25_T2m, "T2m_6h_q25_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q75_T2m, "T2m_6h_q75_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q90_T2m, "T2m_6h_q90_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_min_T2m, "T2m_6h_min_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_max_T2m, "T2m_6h_max_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)

write.table(merge_HRES_DPT, "DPT_6h_HRES_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_avg_DPT, "DPT_6h_mean_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_sd_DPT, "DPT_6h_sd_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_median_DPT, "DPT_6h_median_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q10_DPT, "DPT_6h_q10_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q25_DPT, "DPT_6h_q25_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q75_DPT, "DPT_6h_q75_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q90_DPT, "DPT_6h_q90_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_min_DPT, "DPT_6h_min_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_max_DPT, "DPT_6h_max_2021_Feb.txt", sep = "\t", col.names = T, row.names = F)


Sys.time()










#####（For WSPD only#####################################
#WSPD (m/s)
Sys.time()
setwd('E:/ECMWFdata/ECME_2018')
first_category_name = list.files("Dec")
dir = paste("./Dec/",first_category_name,sep="")
n = length(dir) 
n_sub<-rep(0,n)
n_sub<-as.data.frame(n_sub)
n_sub<-t(n_sub)
head(n_sub)
#create the df for HRES
merge_HRES_WSPD = as.data.frame(matrix(1:40,nrow = 1, ncol = 40))
names(merge_HRES_WSPD) = c(1:40)
merge_HRES_WSPD$second_category <- 'second_category'
merge_HRES_WSPD$first_category <- 'first_category'



#create the df for avg
merge_avg_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_avg_WSPD) = c(1:60)
merge_avg_WSPD$second_category <- 'second_category'
merge_avg_WSPD$first_category <- 'first_category'




#create the df for sd
merge_sd_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_sd_WSPD) = c(1:60)
merge_sd_WSPD$second_category <- 'second_category'
merge_sd_WSPD$first_category <- 'first_category'


#create the df for median
merge_median_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_median_WSPD) = c(1:60)
merge_median_WSPD$second_category <- 'second_category'
merge_median_WSPD$first_category <- 'first_category'



#create the df for q10
merge_q10_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q10_WSPD) = c(1:60)
merge_q10_WSPD$second_category <- 'second_category'
merge_q10_WSPD$first_category <- 'first_category'


#create the df for q25
merge_q25_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q25_WSPD) = c(1:60)
merge_q25_WSPD$second_category <- 'second_category'
merge_q25_WSPD$first_category <- 'first_category'


#createe the df for q75
merge_q75_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q75_WSPD) = c(1:60)
merge_q75_WSPD$second_category <- 'second_category'
merge_q75_WSPD$first_category <- 'first_category'


#create the df for q90
merge_q90_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_q90_WSPD) = c(1:60)
merge_q90_WSPD$second_category <- 'second_category'
merge_q90_WSPD$first_category <- 'first_category'


#create the df for min
merge_min_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_min_WSPD) = c(1:60)
merge_min_WSPD$second_category <- 'second_category'
merge_min_WSPD$first_category <- 'first_category'


#create the df for max
merge_max_WSPD = as.data.frame(matrix(1:60,nrow = 1, ncol = 60))
names(merge_max_WSPD) = c(1:60)
merge_max_WSPD$second_category <- 'second_category'
merge_max_WSPD$first_category <- 'first_category'





for(i in 1:n){         #对于每个一级目录(文件夹)
  b=list.files(dir[i]) #b是列出每个一级目录(文件夹)中每个xlsx文件的名称
  n_sub[i]=length(b)   #得到一级目录(文件夹)下xlsx的文件个数:n_sub
  
  for(j in 1:n_sub[i]){     #对于每个一级目录(文件夹)下的每个xlsx文件
    data<-read.delim(file=paste(dir[i],'/',b[j],sep=''),sep = "", header = F) #读取xlsx文件
    data_WSPD = data[161:213,]
    data_WSPD[,48:83]=trunc(data_WSPD[,48:83]/1000000)
    data_WSPD[2,64]=data_WSPD[2,64]/10^234
    data_WSPD = as.data.frame(lapply(data_WSPD,as.numeric))
    data_WSPD[2:53,]=data_WSPD[2:53,]/10
    WSPD_6hstep = data_WSPD[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]

    
    WSPD_HRES = WSPD_6hstep[2,1:40]

    WSPD_ENS = WSPD_6hstep[3:53,]
    WSPD_ENS = as.data.frame(lapply(WSPD_ENS,as.numeric))

    
    avg_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))

    sd_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))

    
    median_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))

    
    q10_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))

    
    q25_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))

    q75_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))

    
    
    q90_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))

    
    min_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))

    
    max_var_df_WSPD <- as.data.frame(matrix(1:60,nrow = 1, ncol = 60))

    
    
    for(k in 1:60){ 
      avg_var_df_WSPD[,k] <- mean(WSPD_ENS[, k])

      sd_var_df_WSPD[,k] <- sd(WSPD_ENS[, k])

      
      median_var_df_WSPD[,k] <- median(WSPD_ENS[, k])

      q10_var_df_WSPD[,k] <- quantile(WSPD_ENS[,k],probs=0.1)

      
      q25_var_df_WSPD[,k] <- quantile(WSPD_ENS[,k],probs=0.25)

      
      q75_var_df_WSPD[,k] <- quantile(WSPD_ENS[,k],probs=0.75)

      
      q90_var_df_WSPD[,k] <- quantile(WSPD_ENS[,k],probs=0.9)

      
      min_var_df_WSPD[,k] <- min(WSPD_ENS[, k])

      max_var_df_WSPD[,k] <- max(WSPD_ENS[, k])

      
    }    
    #new row for HRES
    names(WSPD_HRES)<-c(1:40)
    WSPD_HRES$second_category<-substr(b[j],14,27)      
    WSPD_HRES$first_category<-first_category_name[i]   
    merge_HRES_WSPD<-rbind(merge_HRES_WSPD,WSPD_HRES)

    #new row for avg
    names(avg_var_df_WSPD)<-c(1:60)
    avg_var_df_WSPD$second_category<-substr(b[j],14,27)       
    avg_var_df_WSPD$first_category<-first_category_name[i] 
    merge_avg_WSPD<-rbind(merge_avg_WSPD,avg_var_df_WSPD)

    
    #new row for sd
    names(sd_var_df_WSPD)<-c(1:60)
    sd_var_df_WSPD$second_category<-substr(b[j],14,27)       
    sd_var_df_WSPD$first_category<-first_category_name[i] 
    merge_sd_WSPD<-rbind(merge_sd_WSPD,sd_var_df_WSPD)
    

    
    #new row for median
    names(median_var_df_WSPD)<-c(1:60)
    median_var_df_WSPD$second_category<-substr(b[j],14,27)       
    median_var_df_WSPD$first_category<-first_category_name[i] 
    merge_median_WSPD<-rbind(merge_median_WSPD,median_var_df_WSPD)
    

    #new row for q10
    names(q10_var_df_WSPD)<-c(1:60)
    q10_var_df_WSPD$second_category<-substr(b[j],14,27)       
    q10_var_df_WSPD$first_category<-first_category_name[i] 
    merge_q10_WSPD<-rbind(merge_q10_WSPD,q10_var_df_WSPD)
    

    
    #new row for q25
    names(q25_var_df_WSPD)<-c(1:60)
    q25_var_df_WSPD$second_category<-substr(b[j],14,27)       
    q25_var_df_WSPD$first_category<-first_category_name[i] 
    merge_q25_WSPD<-rbind(merge_q25_WSPD,q25_var_df_WSPD)
    

    
    #new row for q75
    names(q75_var_df_WSPD)<-c(1:60)
    q75_var_df_WSPD$second_category<-substr(b[j],14,27)       
    q75_var_df_WSPD$first_category<-first_category_name[i] 
    merge_q75_WSPD<-rbind(merge_q75_WSPD,q75_var_df_WSPD)

    
    
    #new row for q90
    names(q90_var_df_WSPD)<-c(1:60)
    q90_var_df_WSPD$second_category<-substr(b[j],14,27)       
    q90_var_df_WSPD$first_category<-first_category_name[i] 
    merge_q90_WSPD<-rbind(merge_q90_WSPD,q90_var_df_WSPD)

    
    #new row for min
    names(min_var_df_WSPD)<-c(1:60)
    min_var_df_WSPD$second_category<-substr(b[j],14,27)       
    min_var_df_WSPD$first_category<-first_category_name[i] 
    merge_min_WSPD<-rbind(merge_min_WSPD,min_var_df_WSPD)
    

    #new row for max
    names(max_var_df_WSPD)<-c(1:60)
    max_var_df_WSPD$second_category<-substr(b[j],14,27)       
    max_var_df_WSPD$first_category<-first_category_name[i] 
    merge_max_WSPD<-rbind(merge_max_WSPD,max_var_df_WSPD)
    

  }
}   

write.table(merge_HRES_WSPD, "WSPD_6h_HRES.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_avg_WSPD, "WSPD_6h_mean.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_sd_WSPD, "WSPD_6h_sd.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_median_WSPD, "WSPD_6h_median.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q10_WSPD, "WSPD_6h_q10.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q25_WSPD, "WSPD_6h_q25.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q75_WSPD, "WSPD_6h_q75.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q90_WSPD, "WSPD_6h_q90.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_min_WSPD, "WSPD_6h_min.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_max_WSPD, "WSPD_6h_max.txt", sep = "\t", col.names = T, row.names = F)




######Temporally match datasets######################################################################
#Read forecast data and radar data


#setwd("E:/ECMWFdata/pre-process/Y2018")
#filelist = list.files(pattern="^RR_6h_q10_2018")
#filelist1 = list.files(pattern = "RR_6h_q10_2018.*")
#setwd("E:/ECMWFdata/pre-process/Y2019")
#filelist2 = list.files(pattern = "RR_6h_q10_2019.*")
#setwd("E:/ECMWFdata/pre-process/Y2020")
#filelist3 = list.files(pattern = "RR_6h_q10_2020.*")
#setwd("E:/ECMWFdata/pre-process/Y2021")
#filelist4 = list.files(pattern = "RR_6h_q10_2021.*")
#file

#data1 = as.data.frame(NULL)

#file1 = paste("E:/ECMWFdata/pre-process/Y2018", filelist1, sep = "/")

#for (i in 1:(length(filelist1))) {
# data = read.table(file1[i],sep = "", header = T)
# data = data[-1,]
# data1 = rbind.data.frame(data1,data)
#}

#file1 = paste("E:/ECMWFdata/pre-process/Y2019", filelist1, sep = "/")

#for (i in 1:(length(filelist1))) {
 # data = read.table(file1[i],sep = "", header = T)
 # data = data[-1,]
 # data1 = rbind.data.frame(data1,data)
#}

#Only select the fourth column of forecast data (for forecasting after 19-24 hours)
library(dplyr)
ECmean201811 = read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_mean_2018_11.txt",sep = "", header = T)
ECmean201811 = ECmean201811[-1,]
ECmean201812 = read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_mean_2018_12.txt",sep = "", header = T)
ECmean201812 = ECmean201812[-1,]
ECmean201901 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_mean_2019_01.txt",sep = "", header = T)
ECmean201901 = ECmean201901[-1,]
ECmean201902 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_mean_2019_02.txt",sep = "", header = T)
ECmean201902 = ECmean201902[-1,]
ECmean201903 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_mean_2019_03.txt",sep = "", header = T)
ECmean201903 = ECmean201903[-1,]
ECmean201904 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_mean_2019_04.txt",sep = "", header = T)
ECmean201904 = ECmean201904[-1,]
ECmean201910 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_mean_2019_10.txt",sep = "", header = T)
ECmean201910 = ECmean201910[-1,]
ECmean201911 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_mean_2019_11.txt",sep = "", header = T)
ECmean201911 = ECmean201911[-1,]
ECmean201912 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_mean_2019_12.txt",sep = "", header = T)
ECmean201912 = ECmean201912[-1,]
ECmean202001 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_mean_2020_01.txt",sep = "", header = T)
ECmean202001 = ECmean202001[-1,]
ECmean202002 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_mean_2020_02.txt",sep = "", header = T)
ECmean202002 = ECmean202002[-1,]
ECmean202003 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_mean_2020_03.txt",sep = "", header = T)
ECmean202003 = ECmean202003[-1,]
ECmean202004 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_mean_2020_04.txt",sep = "", header = T)
ECmean202004 = ECmean202004[-1,]
ECmean202010 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_mean_2020_10.txt",sep = "", header = T)
ECmean202010 = ECmean202010[-1,]
ECmean202011 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_mean_2020_11.txt",sep = "", header = T)
ECmean202011 = ECmean202011[-1,]
ECmean202012 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_mean_2020_12.txt",sep = "", header = T)
ECmean202012 = ECmean202012[-1,]
ECmean202101 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_mean_2021_01.txt",sep = "", header = T)
ECmean202101 = ECmean202101[-1,]
ECmean202102 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_mean_2021_02.txt",sep = "", header = T)
ECmean202102 = ECmean202102[-1,]
ECmean202103 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_mean_2021_03.txt",sep = "", header = T)
ECmean202103 = ECmean202103[-1,]
ECmean202104 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_mean_2021_04.txt",sep = "", header = T)
ECmean202104 = ECmean202104[-1,]

ECmean = rbind.data.frame(ECmean201811, ECmean201812, ECmean201901, ECmean201902, ECmean201903, ECmean201904, ECmean201910,
                          ECmean201911, ECmean201912, ECmean202001, ECmean202002, ECmean202003, ECmean202004, ECmean202010, 
                          ECmean202011, ECmean202012, ECmean202101, ECmean202102, ECmean202103, ECmean202104)







ECsd201811 = read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_sd_2018_11.txt",sep = "", header = T)
ECsd201811 = ECsd201811[-1,]
ECsd201812 = read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_sd_2018_12.txt",sep = "", header = T)
ECsd201812 = ECsd201812[-1,]
ECsd201901 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_sd_2019_01.txt",sep = "", header = T)
ECsd201901 = ECsd201901[-1,]
ECsd201902 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_sd_2019_02.txt",sep = "", header = T)
ECsd201902 = ECsd201902[-1,]
ECsd201903 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_sd_2019_03.txt",sep = "", header = T)
ECsd201903 = ECsd201903[-1,]
ECsd201904 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_sd_2019_04.txt",sep = "", header = T)
ECsd201904 = ECsd201904[-1,]
ECsd201910 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_sd_2019_10.txt",sep = "", header = T)
ECsd201910 = ECsd201910[-1,]
ECsd201911 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_sd_2019_11.txt",sep = "", header = T)
ECsd201911 = ECsd201911[-1,]
ECsd201912 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_sd_2019_12.txt",sep = "", header = T)
ECsd201912 = ECsd201912[-1,]
ECsd202001 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_sd_2020_01.txt",sep = "", header = T)
ECsd202001 = ECsd202001[-1,]
ECsd202002 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_sd_2020_02.txt",sep = "", header = T)
ECsd202002 = ECsd202002[-1,]
ECsd202003 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_sd_2020_03.txt",sep = "", header = T)
ECsd202003 = ECsd202003[-1,]
ECsd202004 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_sd_2020_04.txt",sep = "", header = T)
ECsd202004 = ECsd202004[-1,]
ECsd202010 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_sd_2020_10.txt",sep = "", header = T)
ECsd202010 = ECsd202010[-1,]
ECsd202011 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_sd_2020_11.txt",sep = "", header = T)
ECsd202011 = ECsd202011[-1,]
ECsd202012 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_sd_2020_12.txt",sep = "", header = T)
ECsd202012 = ECsd202012[-1,]
ECsd202101 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_sd_2021_01.txt",sep = "", header = T)
ECsd202101 = ECsd202101[-1,]
ECsd202102 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_sd_2021_02.txt",sep = "", header = T)
ECsd202102 = ECsd202102[-1,]
ECsd202103 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_sd_2021_03.txt",sep = "", header = T)
ECsd202103 = ECsd202103[-1,]
ECsd202104 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_sd_2021_04.txt",sep = "", header = T)
ECsd202104 = ECsd202104[-1,]

ECsd = rbind.data.frame(ECsd201811, ECsd201812, ECsd201901, ECsd201902, ECsd201903, ECsd201904, ECsd201910, ECsd201911, ECsd201912,
                        ECsd202001, ECsd202002, ECsd202003, ECsd202004, ECsd202010, ECsd202011, ECsd202012, ECsd202101, ECsd202102,
                        ECsd202103, ECsd202104)


ECmedian201811 = read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_median_2018_11.txt",sep = "", header = T)
ECmedian201811 = ECmedian201811[-1,]
ECmedian201812 = read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_median_2018_12.txt",sep = "", header = T)
ECmedian201812 = ECmedian201812[-1,]
ECmedian201901 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_median_2019_01.txt",sep = "", header = T)
ECmedian201901 = ECmedian201901[-1,]
ECmedian201902 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_median_2019_02.txt",sep = "", header = T)
ECmedian201902 = ECmedian201902[-1,]
ECmedian201903 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_median_2019_03.txt",sep = "", header = T)
ECmedian201903 = ECmedian201903[-1,]
ECmedian201904 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_median_2019_04.txt",sep = "", header = T)
ECmedian201904 = ECmedian201904[-1,]
ECmedian201910 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_median_2019_10.txt",sep = "", header = T)
ECmedian201910 = ECmedian201910[-1,]
ECmedian201911 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_median_2019_11.txt",sep = "", header = T)
ECmedian201911 = ECmedian201911[-1,]
ECmedian201912 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_median_2019_12.txt",sep = "", header = T)
ECmedian201912 = ECmedian201912[-1,]
ECmedian202001 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_median_2020_01.txt",sep = "", header = T)
ECmedian202001 = ECmedian202001[-1,]
ECmedian202002 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_median_2020_02.txt",sep = "", header = T)
ECmedian202002 = ECmedian202002[-1,]
ECmedian202003 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_median_2020_03.txt",sep = "", header = T)
ECmedian202003 = ECmedian202003[-1,]
ECmedian202004 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_median_2020_04.txt",sep = "", header = T)
ECmedian202004 = ECmedian202004[-1,]
ECmedian202010 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_median_2020_10.txt",sep = "", header = T)
ECmedian202010 = ECmedian202010[-1,]
ECmedian202011 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_median_2020_11.txt",sep = "", header = T)
ECmedian202011 = ECmedian202011[-1,]
ECmedian202012 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_median_2020_12.txt",sep = "", header = T)
ECmedian202012 = ECmedian202012[-1,]
ECmedian202101 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_median_2021_01.txt",sep = "", header = T)
ECmedian202101 = ECmedian202101[-1,]
ECmedian202102 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_median_2021_02.txt",sep = "", header = T)
ECmedian202102 = ECmedian202102[-1,]
ECmedian202103 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_median_2021_03.txt",sep = "", header = T)
ECmedian202103 = ECmedian202103[-1,]
ECmedian202104 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_median_2021_04.txt",sep = "", header = T)
ECmedian202104 = ECmedian202104[-1,]

ECmedian = rbind.data.frame(ECmedian201811,ECmedian201812, ECmedian201901, ECmedian201902, ECmedian201903, ECmedian201904, ECmedian201910, ECmedian201911,
                            ECmedian201912, ECmedian202001, ECmedian202002, ECmedian202003, ECmedian202004, ECmedian202010, ECmedian202011, ECmedian202012,
                            ECmedian202101, ECmedian202102, ECmedian202103, ECmedian202104)

#setwd("E:/ECMWFdata/pre-process/Y2018")
#filelist = list.files(pattern="^RR_6h_q10_2018")
#filelist1 = list.files(pattern = "RR_6h_q10_2018.*")
#setwd("E:/ECMWFdata/pre-process/Y2019")
#filelist2 = list.files(pattern = "RR_6h_q10_2019.*")
#setwd("E:/ECMWFdata/pre-process/Y2020")
#filelist3 = list.files(pattern = "RR_6h_q10_2020.*")
#setwd("E:/ECMWFdata/pre-process/Y2021")
#filelist4 = list.files(pattern = "RR_6h_q10_2021.*")
#file

#data1 = as.data.frame(NULL)

#file1 = paste("E:/ECMWFdata/pre-process/Y2018", filelist1, sep = "/")

#for (i in 1:(length(filelist1))) {
# data = read.table(file1[i],sep = "", header = T)
# data = data[-1,]
# data1 = rbind.data.frame(data1,data)
#}

#file1 = paste("E:/ECMWFdata/pre-process/Y2019", filelist1, sep = "/")

#for (i in 1:(length(filelist1))) {
# data = read.table(file1[i],sep = "", header = T)
# data = data[-1,]
# data1 = rbind.data.frame(data1,data)
#}







ECq10_201811 = read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_q10_2018_11.txt",sep = "", header = T)
ECq10_201811 = ECq10_201811[-1,]
ECq10_201812 = read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_q10_2018_12.txt",sep = "", header = T)
ECq10_201812 = ECq10_201812[-1,]
ECq10_201901 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q10_2019_01.txt",sep = "", header = T)
ECq10_201901 = ECq10_201901[-1,]
ECq10_201902 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q10_2019_02.txt",sep = "", header = T)
ECq10_201902 = ECq10_201902[-1,]
ECq10_201903 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q10_2019_03.txt",sep = "", header = T)
ECq10_201903 = ECq10_201903[-1,]
ECq10_201904 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q10_2019_04.txt",sep = "", header = T)
ECq10_201904 = ECq10_201904[-1,]
ECq10_201910 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q10_2019_10.txt",sep = "", header = T)
ECq10_201910 = ECq10_201910[-1,]
ECq10_201911 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q10_2019_11.txt",sep = "", header = T)
ECq10_201911 = ECq10_201911[-1,]
ECq10_201912 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q10_2019_12.txt",sep = "", header = T)
ECq10_201912 = ECq10_201912[-1,]
ECq10_202001 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q10_2020_01.txt",sep = "", header = T)
ECq10_202001 = ECq10_202001[-1,]
ECq10_202002 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q10_2020_02.txt",sep = "", header = T)
ECq10_202002 = ECq10_202002[-1,]
ECq10_202003 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q10_2020_03.txt",sep = "", header = T)
ECq10_202003 = ECq10_202003[-1,]
ECq10_202004 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q10_2020_04.txt",sep = "", header = T)
ECq10_202004 = ECq10_202004[-1,]
ECq10_202010 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q10_2020_10.txt",sep = "", header = T)
ECq10_202010 = ECq10_202010[-1,]
ECq10_202011 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q10_2020_11.txt",sep = "", header = T)
ECq10_202011 = ECq10_202011[-1,]
ECq10_202012 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q10_2020_12.txt",sep = "", header = T)
ECq10_202012 = ECq10_202012[-1,]
ECq10_202101 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q10_2021_01.txt",sep = "", header = T)
ECq10_202101 = ECq10_202101[-1,]
ECq10_202102 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q10_2021_02.txt",sep = "", header = T)
ECq10_202102 = ECq10_202102[-1,]
ECq10_202103 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q10_2021_03.txt",sep = "", header = T)
ECq10_202103 = ECq10_202103[-1,]
ECq10_202104 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q10_2021_04.txt",sep = "", header = T)
ECq10_202104 = ECq10_202104[-1,]


ECq90_201811 = read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_q90_2018_11.txt",sep = "", header = T)
ECq90_201811 = ECq90_201811[-1,]
ECq90_201812 = read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_q90_2018_12.txt",sep = "", header = T)
ECq90_201812 = ECq90_201812[-1,]
ECq90_201901 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q90_2019_01.txt",sep = "", header = T)
ECq90_201901 = ECq90_201901[-1,]
ECq90_201902 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q90_2019_02.txt",sep = "", header = T)
ECq90_201902 = ECq90_201902[-1,]
ECq90_201903 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q90_2019_03.txt",sep = "", header = T)
ECq90_201903 = ECq90_201903[-1,]
ECq90_201904 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q90_2019_04.txt",sep = "", header = T)
ECq90_201904 = ECq90_201904[-1,]
ECq90_201910 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q90_2019_10.txt",sep = "", header = T)
ECq90_201910 = ECq90_201910[-1,]
ECq90_201911 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q90_2019_11.txt",sep = "", header = T)
ECq90_201911 = ECq90_201911[-1,]
ECq90_201912 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_q90_2019_12.txt",sep = "", header = T)
ECq90_201912 = ECq90_201912[-1,]
ECq90_202001 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q90_2020_01.txt",sep = "", header = T)
ECq90_202001 = ECq90_202001[-1,]
ECq90_202002 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q90_2020_02.txt",sep = "", header = T)
ECq90_202002 = ECq90_202002[-1,]
ECq90_202003 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q90_2020_03.txt",sep = "", header = T)
ECq90_202003 = ECq90_202003[-1,]
ECq90_202004 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q90_2020_04.txt",sep = "", header = T)
ECq90_202004 = ECq90_202004[-1,]
ECq90_202010 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q90_2020_10.txt",sep = "", header = T)
ECq90_202010 = ECq90_202010[-1,]
ECq90_202011 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q90_2020_11.txt",sep = "", header = T)
ECq90_202011 = ECq90_202011[-1,]
ECq90_202012 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_q90_2020_12.txt",sep = "", header = T)
ECq90_202012 = ECq90_202012[-1,]
ECq90_202101 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q90_2021_01.txt",sep = "", header = T)
ECq90_202101 = ECq90_202101[-1,]
ECq90_202102 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q90_2021_02.txt",sep = "", header = T)
ECq90_202102 = ECq90_202102[-1,]
ECq90_202103 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q90_2021_03.txt",sep = "", header = T)
ECq90_202103 = ECq90_202103[-1,]
ECq90_202104 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_q90_2021_04.txt",sep = "", header = T)
ECq90_202104 = ECq90_202104[-1,]




ECp0_201811 = read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_p0_2018_11.txt",sep = "", header = T)
ECp0_201811 = ECp0_201811[-1,]
ECp0_201812 = read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_p0_2018_12.txt",sep = "", header = T)
ECp0_201812 = ECp0_201812[-1,]
ECp0_201901 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p0_2019_01.txt",sep = "", header = T)
ECp0_201901 = ECp0_201901[-1,]
ECp0_201902 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p0_2019_02.txt",sep = "", header = T)
ECp0_201902 = ECp0_201902[-1,]
ECp0_201903 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p0_2019_03.txt",sep = "", header = T)
ECp0_201903 = ECp0_201903[-1,]
ECp0_201904 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p0_2019_04.txt",sep = "", header = T)
ECp0_201904 = ECp0_201904[-1,]
ECp0_201910 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p0_2019_10.txt",sep = "", header = T)
ECp0_201910 = ECp0_201910[-1,]
ECp0_201911 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p0_2019_11.txt",sep = "", header = T)
ECp0_201911 = ECp0_201911[-1,]
ECp0_201912 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_p0_2019_12.txt",sep = "", header = T)
ECp0_201912 = ECp0_201912[-1,]
ECp0_202001 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_p0_2020_01.txt",sep = "", header = T)
ECp0_202001 = ECp0_202001[-1,]
ECp0_202002 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_p0_2020_02.txt",sep = "", header = T)
ECp0_202002 = ECp0_202002[-1,]
ECp0_202003 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_p0_2020_03.txt",sep = "", header = T)
ECp0_202003 = ECp0_202003[-1,]
ECp0_202004 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_p0_2020_04.txt",sep = "", header = T)
ECp0_202004 = ECp0_202004[-1,]
ECp0_202010 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_p0_2020_10.txt",sep = "", header = T)
ECp0_202010 = ECp0_202010[-1,]
ECp0_202011 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_p0_2020_11.txt",sep = "", header = T)
ECp0_202011 = ECp0_202011[-1,]
ECp0_202012 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_p0_2020_12.txt",sep = "", header = T)
ECp0_202012 = ECp0_202012[-1,]
ECp0_202101 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_p0_2021_01.txt",sep = "", header = T)
ECp0_202101 = ECp0_202101[-1,]
ECp0_202102 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_p0_2021_02.txt",sep = "", header = T)
ECp0_202102 = ECp0_202102[-1,]
ECp0_202103 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_p0_2021_03.txt",sep = "", header = T)
ECp0_202103 = ECp0_202103[-1,]
ECp0_202104 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_p0_2021_04.txt",sep = "", header = T)
ECp0_202104 = ECp0_202104[-1,]




ECq10 = rbind.data.frame(ECq10_201811,ECq10_201812, ECq10_201901, ECq10_201902, ECq10_201903, ECq10_201904, ECq10_201910, ECq10_201911,
                         ECq10_201912, ECq10_202001, ECq10_202002, ECq10_202003, ECq10_202004, ECq10_202010, ECq10_202011, ECq10_202012,
                         ECq10_202101, ECq10_202102, ECq10_202103, ECq10_202104)
#ECq10 = ECq10[-(66551:73810),]
#ECq10 = ECq10[-(159237:166496),]


ECq90 = rbind.data.frame(ECq90_201811,ECq90_201812, ECq90_201901, ECq90_201902, ECq90_201903, ECq90_201904, ECq90_201910, ECq90_201911,
                         ECq90_201912, ECq90_202001, ECq90_202002, ECq90_202003, ECq90_202004, ECq90_202010, ECq90_202011, ECq90_202012,
                         ECq90_202101, ECq90_202102, ECq90_202103, ECq90_202104)

#ECq90 = ECq90[-(66551:73810),]
#ECq90 = ECq90[-(159237:166496),]


ECp0 = rbind.data.frame(ECp0_201811,ECp0_201812, ECp0_201901, ECp0_201902, ECp0_201903, ECp0_201904, ECp0_201910, ECp0_201911,
                        ECp0_201912, ECp0_202001, ECp0_202002, ECp0_202003, ECp0_202004, ECp0_202010, ECp0_202011, ECp0_202012,
                        ECp0_202101, ECp0_202102, ECp0_202103, ECp0_202104)
#ECp0 = ECp0[-(66551:73810),]
#ECp0 = ECp0[-(159237:166496),]


rm(ECmean201811,ECmean201812,ECmean201901,ECmean201902,ECmean201903,ECmean201904,ECmean201910,ECmean201911,
   ECmean201912,ECmean202001,ECmean202002,ECmean202003,ECmean202004,ECmean202010,ECmean202011,ECmean202012,
   ECmean202101,ECmean202102,ECmean202103,ECmean202104,ECsd201811,ECsd201812,ECsd201901,ECsd201902,ECsd201903,
   ECsd201904,ECsd201910,ECsd201911,ECsd201912,ECsd202001,ECsd202002,ECsd202003,ECsd202004,ECsd202010,
   ECsd202011,ECsd202012,ECsd202101,ECsd202102,ECsd202103,ECsd202104)  
rm(ECmedian201811,ECmedian201812, ECmedian201901, ECmedian201902, ECmedian201903, ECmedian201904, ECmedian201910, ECmedian201911,
   ECmedian201912, ECmedian202001, ECmedian202002, ECmedian202003, ECmedian202004, ECmedian202010, ECmedian202011, ECmedian202012,
   ECmedian202101, ECmedian202102, ECmedian202103, ECmedian202104,ECq10_201811,ECq10_201812, ECq10_201901, ECq10_201902, ECq10_201903, ECq10_201904, ECq10_201910, ECq10_201911,
   ECq10_201912, ECq10_202001, ECq10_202002, ECq10_202003, ECq10_202004, ECq10_202010, ECq10_202011, ECq10_202012,
   ECq10_202101, ECq10_202102, ECq10_202103, ECq10_202104,ECq90_201811,ECq90_201812, ECq90_201901, ECq90_201902, ECq90_201903, ECq90_201904, ECq90_201910, ECq90_201911,
   ECq90_201912, ECq90_202001, ECq90_202002, ECq90_202003, ECq90_202004, ECq90_202010, ECq90_202011, ECq90_202012,
   ECq90_202101, ECq90_202102, ECq90_202103, ECq90_202104,ECp0_201811,ECp0_201812, ECp0_201901, ECp0_201902, ECp0_201903, ECp0_201904, ECp0_201910, ECp0_201911,
   ECp0_201912, ECp0_202001, ECp0_202002, ECp0_202003, ECp0_202004, ECp0_202010, ECp0_202011, ECp0_202012,
   ECp0_202101, ECp0_202102, ECp0_202103, ECp0_202104)



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







ECmean$grid = as.numeric(substr(ECmean$second_category,12,15))
ECsd$grid = as.numeric(substr(ECsd$second_category,12,15))
ECmedian$grid = as.numeric(substr(ECmedian$second_category,12,15))
ECq10$grid = as.numeric(substr(ECq10$second_category,12,15))
ECq90$grid = as.numeric(substr(ECq90$second_category,12,15))
ECp0$grid = as.numeric(substr(ECp0$second_category,12,15))



ECmean$dftime = as.POSIXct("2000-1-1 01:00:00", tz = "UTC")
year(ECmean$dftime) = as.numeric(substr(ECmean$first_category,10,13))
month(ECmean$dftime) = as.numeric(substr(ECmean$first_category,14,15))
mday(ECmean$dftime) = as.numeric(substr(ECmean$first_category,16,17))
hour(ECmean$dftime) = as.numeric(substr(ECmean$first_category,18,19))
minute(ECmean$dftime) = 00
second(ECmean$dftime) = 00
ECmean$dftime = as.character(ECmean$dftime)



ECsd$dftime = as.POSIXct("2000-1-1 01:00:00", tz = "UTC")
year(ECsd$dftime) = as.numeric(substr(ECsd$first_category,10,13))
month(ECsd$dftime) = as.numeric(substr(ECsd$first_category,14,15))
mday(ECsd$dftime) = as.numeric(substr(ECsd$first_category,16,17))
hour(ECsd$dftime) = as.numeric(substr(ECsd$first_category,18,19))
minute(ECsd$dftime) = 00
second(ECsd$dftime) = 00
ECsd$dftime = as.character(ECsd$dftime)


ECmedian$dftime = as.POSIXct("2000-1-1 01:00:00", tz = "UTC")
year(ECmedian$dftime) = as.numeric(substr(ECmedian$first_category,10,13))
month(ECmedian$dftime) = as.numeric(substr(ECmedian$first_category,14,15))
mday(ECmedian$dftime) = as.numeric(substr(ECmedian$first_category,16,17))
hour(ECmedian$dftime) = as.numeric(substr(ECmedian$first_category,18,19))
minute(ECmedian$dftime) = 00
second(ECmedian$dftime) = 00
ECmedian$dftime = as.character(ECmedian$dftime)


ECq10$dftime = as.POSIXct("2000-1-1 01:00:00", tz = "UTC")
year(ECq10$dftime) = as.numeric(substr(ECq10$first_category,10,13))
month(ECq10$dftime) = as.numeric(substr(ECq10$first_category,14,15))
mday(ECq10$dftime) = as.numeric(substr(ECq10$first_category,16,17))
hour(ECq10$dftime) = as.numeric(substr(ECq10$first_category,18,19))
minute(ECq10$dftime) = 00
second(ECq10$dftime) = 00
ECq10$dftime = as.character(ECq10$dftime)

ECq90$dftime = as.POSIXct("2000-1-1 01:00:00", tz = "UTC")
year(ECq90$dftime) = as.numeric(substr(ECq90$first_category,10,13))
month(ECq90$dftime) = as.numeric(substr(ECq90$first_category,14,15))
mday(ECq90$dftime) = as.numeric(substr(ECq90$first_category,16,17))
hour(ECq90$dftime) = as.numeric(substr(ECq90$first_category,18,19))
minute(ECq90$dftime) = 00
second(ECq90$dftime) = 00
ECq90$dftime = as.character(ECq90$dftime)

ECp0$dftime = as.POSIXct("2000-1-1 01:00:00", tz = "UTC")
year(ECp0$dftime) = as.numeric(substr(ECp0$first_category,10,13))
month(ECp0$dftime) = as.numeric(substr(ECp0$first_category,14,15))
mday(ECp0$dftime) = as.numeric(substr(ECp0$first_category,16,17))
hour(ECp0$dftime) = as.numeric(substr(ECp0$first_category,18,19))
minute(ECp0$dftime) = 00
second(ECp0$dftime) = 00
ECp0$dftime = as.character(ECp0$dftime)


ECsub = ECmean[,61:64]
Rsub = Rdata[,1:2]
ECmatch = left_join(ECsub,Rsub,by = c("dftime"="Rtime"))
ECmean = cbind.data.frame(ECmean,ECmatch$RmatchID)
ECmean = ECmean[which(ECmean$`ECmatch$RmatchID` <665 | ECmean$`ECmatch$RmatchID` >783 & ECmean$`ECmatch$RmatchID`<1517
                      | ECmean$`ECmatch$RmatchID`>1635 & ECmean$`ECmatch$RmatchID`<2365 ),]
ECsd = cbind.data.frame(ECsd,ECmatch$RmatchID)
ECsd = ECsd[which(ECsd$`ECmatch$RmatchID` <665 | ECsd$`ECmatch$RmatchID` >783  & ECsd$`ECmatch$RmatchID`<1517
                  | ECsd$`ECmatch$RmatchID`>1635 & ECsd$`ECmatch$RmatchID`<2365 ),]
ECmedian = cbind.data.frame(ECmedian,ECmatch$RmatchID)
ECmedian = ECmedian[which(ECmedian$`ECmatch$RmatchID` <665 | ECmedian$`ECmatch$RmatchID` >783 & ECmedian$`ECmatch$RmatchID`<1517
                          | ECmedian$`ECmatch$RmatchID`>1635 & ECmedian$`ECmatch$RmatchID`<2365 ),]
ECq10 = cbind.data.frame(ECq10,ECmatch$RmatchID)
ECq10 = ECq10[which(ECq10$`ECmatch$RmatchID` <665 | ECq10$`ECmatch$RmatchID` >783  & ECq10$`ECmatch$RmatchID`<1517
                    | ECq10$`ECmatch$RmatchID`>1635 & ECq10$`ECmatch$RmatchID`<2365 ),]
ECq90 = cbind.data.frame(ECq90,ECmatch$RmatchID)
ECq90 = ECq90[which(ECq90$`ECmatch$RmatchID` <665 | ECq90$`ECmatch$RmatchID` >783  & ECq90$`ECmatch$RmatchID`<1517
                    | ECq90$`ECmatch$RmatchID`>1635 & ECq90$`ECmatch$RmatchID`<2365 ),]
ECp0 = cbind.data.frame(ECp0,ECmatch$RmatchID)
ECp0 = ECp0[which(ECp0$`ECmatch$RmatchID` <665 | ECp0$`ECmatch$RmatchID` >783  & ECp0$`ECmatch$RmatchID`<1517
                  | ECp0$`ECmatch$RmatchID`>1635 & ECp0$`ECmatch$RmatchID`<2365 ),]




#Take an example for predicting （6-12）hours:
inputmean = as.data.frame(ECmean[,c(2,61:65)])
inputsd = as.data.frame(ECsd[,c(2,61:65)])
inputmedian = as.data.frame(ECmedian[,c(2,61:65)])
inputq10 = as.data.frame(ECq10[,c(2,61:65)])
inputq90 = as.data.frame(ECq90[,c(2,61:65)])
inputp0 = as.data.frame(ECp0[,c(2,61:65)])



#  Rdata1 = Rdata[,-c(2425,2426)]
#   Rdata1 = Rdata1[,-c(1:3)]
#  upscal_match = read.table("E:/KNMI_DATA/upscal_match_9km.txt", sep = "", header = T)

EC4meanfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
EC4meanfinal <- EC4meanfinal[-(1:2),]
EC4time <- as.data.frame(matrix(upscal_match$rID,nrow=601,ncol=1))
EC4time <- EC4time[-(1:2),]
EC4sdfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
EC4sdfinal <- EC4sdfinal[-(1:2),]
EC4medianfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
EC4medianfinal <- EC4medianfinal[-(1:2),]
EC4q10final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
EC4q10final <- EC4q10final[-(1:2),]
EC4q90final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
EC4q90final <- EC4q90final[-(1:2),]
EC4p0final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
EC4p0final <- EC4p0final[-(1:2),]
R4final <- as.data.frame(matrix(upscal_match$rID,nrow = 601, ncol = 1))
R4final <- R4final[-(1:2),]
EC4ID <- as.data.frame(matrix(upscal_match$rID,nrow = 601,ncol=1))
EC4ID <- EC4ID[-(1:2),]
R4ID <- as.data.frame(matrix(upscal_match$rID,nrow = 601,ncol=1))
R4ID <- R4ID[-(1:2),]

# Rdata1 = Rdata
# time = as.data.frame(Rdata1$time)
# time = as.data.frame(time[-(1:2),])
# colnames(time) = "time"
#time$time = ymd_hms(time$time)

#   Rdata$Rtime = ymd_hms(Rdata$Rtime)



#  test1 = as.POSIXct(test1, origin = "1970-1-1 00:00:00", tz = "UTC")

#make df for training  

#for (i in unique(inputmean$dftime)){


test = unique(ECmean$`ECmatch$RmatchID`)

for (i in test){
  #  subdata = Rdata[which(Rdata$Rtime == i+hours(24)),]   
  subdata = Rdata[which(Rdata$RmatchID == i+1),]
  
  #    rID = upscal_match[,1]
  #    dfID = upscal_match[,2]   
  #    subdata = rbind.data.frame(rID,dfID,subdata)   
  subdata = t(subdata)
  subdata = as.data.frame(subdata)
  colnames(subdata) = c(subdata[2,1])
  submean = inputmean[which(inputmean$`ECmatch$RmatchID` == i),]
  subsd = inputsd[which(inputsd$`ECmatch$RmatchID` == i),]
  submedian = inputmedian[which(inputmedian$`ECmatch$RmatchID` == i),]
  subq10 = inputq10[which(inputq10$`ECmatch$RmatchID` == i),]
  subq90 = inputq90[which(inputq90$`ECmatch$RmatchID` == i),]
  subp0 = inputp0[which(inputp0$`ECmatch$RmatchID` == i),]
  
  
  
  
  subpre = cbind.data.frame(submean[,1],subsd[,1],submedian[,1],subq10[,1],subq90[,1],subp0)
  colnames(subpre)[1] = "submean"
  colnames(subpre)[2] = "subsd"
  colnames(subpre)[3] = "submedian"
  colnames(subpre)[4] = "subq10"
  colnames(subpre)[5] = "subq90"
  colnames(subpre)[6] = "subp0"
  
  subpre = subpre[,c(1,2,3,4,5,6,9,10)]
  subdata$dfID=as.numeric(upscal_match$dfID)
  subdata$rID=as.numeric(upscal_match$rID)
  match = left_join(subdata,subpre,by = c("dfID"="grid"))
  match = match[-(1:2),]
  EC4Rmean = match[,4]
  EC4sd = match[,5]
  EC4median = match[,6]
  EC4q10 = match[,7]
  EC4q90 = match[,8]
  EC4p0 = match[,9]
  R4match = match[,1]
  dftime = match[,10]
  ECID = match[,2]
  RID = match[,3]
  EC4meanfinal=cbind.data.frame(EC4meanfinal,EC4Rmean)
  EC4sdfinal = cbind.data.frame(EC4sdfinal,EC4sd)
  EC4medianfinal=cbind.data.frame(EC4medianfinal,EC4median)
  EC4q10final=cbind.data.frame(EC4q10final,EC4q10)
  EC4q90final=cbind.data.frame(EC4q90final,EC4q90)
  EC4p0final=cbind.data.frame(EC4p0final,EC4p0)
  R4final = cbind.data.frame(R4final,R4match)
  EC4time <- cbind.data.frame(EC4time,dftime)
  EC4ID = cbind.data.frame(EC4ID,ECID)
  R4ID = cbind.data.frame(R4ID,RID)
}



EC4meanfinal1 = EC4meanfinal[,-1]
EC4meanfinal1 = as.matrix(EC4meanfinal1)
dim(EC4meanfinal1) <- c(599*955,1)
EC4sdfinal1 = EC4sdfinal[,-1]
EC4sdfinal1 = as.matrix(EC4sdfinal1)
dim(EC4sdfinal1) <- c(599*955,1)
EC4medianfinal1 = EC4medianfinal[,-1]
EC4medianfinal1 = as.matrix(EC4medianfinal1)
dim(EC4medianfinal1) <- c(599*955,1)
EC4q10final1 = EC4q10final[,-1]
EC4q10final1 = as.matrix(EC4q10final1)
dim(EC4q10final1) <- c(599*955,1)
EC4q90final1 = EC4q90final[,-1]
EC4q90final1 = as.matrix(EC4q90final1)
dim(EC4q90final1) <- c(599*955,1)
EC4p0final1 = EC4p0final[,-1]
EC4p0final1 = as.matrix(EC4p0final1)
dim(EC4p0final1) <- c(599*955,1)
R4final1 = R4final[,-1]
R4final1 = as.matrix(R4final1)
dim(R4final1) <- c(599*955,1)
EC4time1 = EC4time[,-1]
EC4time1 = as.matrix(EC4time1)
dim(EC4time1) <- c(599*955,1)
EC4ID1 = EC4ID[,-1]
EC4ID1 = as.matrix(EC4ID1)
dim(EC4ID1) <- c(599*955,1)
R4ID1 = R4ID[,-1]
R4ID1 = as.matrix(R4ID1)
dim(R4ID1) <- c(599*955,1)


final = cbind.data.frame(EC4meanfinal1,EC4sdfinal1,EC4medianfinal1,EC4q10final1,EC4q90final1,EC4p0final1,R4final1,EC4time1,EC4ID1,R4ID1)
final[which(final$EC4meanfinal1 <0 ),1]<- 0


rm(dfID,EC4ID,EC4ID1,EC4meanfinal,EC4meanfinal1,EC4sdfinal,EC4sdfinal1,EC4time,EC4time1,
   ECmean,ECsd,inputmean,inputsd,match,new,R4final,R4final1,Rdata,rID,subdata,submean,subpre,subsd,upscal_match)

rm(EC4medianfinal,EC4medianfinal1,EC4p0final,EC4p0final1,EC4q10final,EC4q10final1,EC4q90final,EC4q90final1,
   ECmedian,ECp0,ECq10,ECq90,inputmedian,inputp0,inputq10,inputq90,submedian,subp0,subq10,subq90)


final$EC4meanfinal1 = as.numeric(final$EC4meanfinal1)
final$EC4sdfinal1 = as.numeric(final$EC4sdfinal1)
final$EC4medianfinal1 = as.numeric(final$EC4medianfinal1)
final$EC4q10final1 = as.numeric(final$EC4q10final1)
final$EC4q90final1 = as.numeric(final$EC4q90final1)
final$EC4p0final1 = as.numeric(final$EC4p0final1)

final$R4final1 = as.numeric(final$R4final1)
final = na.omit(final)




library(gamlss)
curt <- function(value) return(value^(1.0/3))






mod0<-gamlss(R4final1~1, data=final, family=ZAGA)

mod1 = stepGAICAll.A(mod0, scope=list(lower=~1,upper=~EC4meanfinal1+EC4sdfinal1+EC4medianfinal1+EC4q10final1+EC4q90final+EC4p0final1), 
                     sigma.scope=list(lower=~1,upper=~EC4meanfinal1+EC4sdfinal1+EC4medianfinal1+EC4q10final1+EC4q90final+EC4p0final1), 
                     nu.scope=list(lower=~1,upper=~EC4meanfinal1+EC4sdfinal1+EC4medianfinal1+EC4q10final1+EC4q90final+EC4p0final1),
                     #direction = "forward",
                     steps = 2) 

saveRDS(object = mod1, file = "toy model for RR.rds")
mod1 = readRDS(file = "toy model for RR.rds")



prd  <- predictAll(mod1,newdata=final[1:20,])
prob <- 1:51 / 52
new_data = as.data.frame(matrix(NA,nrow = 51, ncol = 20))
for (i in 1:20){
new_data[,i]=sapply(prob, function(u) qZAGA(u, mu=prd$mu[i], sigma=prd$sigma[i], nu=prd$nu[i]))
}




CRPS; 


qts  <- switch(env$prm, T2m=sapply(prob, function(u) qNO(u, mu=prd$mu, sigma=prd$sigma)),
               S10m=sapply(prob, function(u) qBCT(u, mu=prd$mu, sigma=prd$sigma, nu=prd$nu, tau=prd$tau)),
               AccPcp1h=sapply(prob, function(u) qZAGA(u, mu=prd$mu, sigma=prd$sigma, nu=prd$nu)))
cat("Range of predicted quantiles:", range(qts), "\n")


mod1 = stepGAICAll.A(mod0, scope=list(lower=~1,upper=~EC4meanfinal1+EC4sdfinal1+EC4medianfinal1+EC4q10final1+EC4q90final1+EC4p0final1), 
                     sigma.scope=list(lower=~1,upper=~EC4meanfinal1+EC4sdfinal1+EC4medianfinal1+EC4q10final1+EC4q90final1+EC4p0final1), 
                     nu.scope=list(lower=~1,upper=~EC4meanfinal1+EC4sdfinal1+EC4medianfinal1+EC4q10final1+EC4q90final1+EC4p0final1)) 


#Mu Coefficients:
#  (Intercept)   EC4q90final1    EC4p0final1    EC4sdfinal1  EC4meanfinal1   EC4q10final1  
#-3.30495        0.14569        2.85005        0.08243        0.23427       -0.08749  
#Sigma Coefficients:
#  (Intercept)     EC4q90final1      EC4p0final1     EC4q10final1      EC4sdfinal1  EC4medianfinal1    EC4meanfinal1  
#0.58493         -0.08516         -0.17233         -0.01741          0.10119          0.08296         -0.11659  
#Nu Coefficients:
#  (Intercept)     EC4q90final1      EC4p0final1      EC4sdfinal1     EC4q10final1  EC4medianfinal1    EC4meanfinal1  
#2.4123          -0.4389          -2.0214          -3.7861          -1.8932           0.6710           0.4136 

qZAGA(0.975, mu = -3.30495+0.14569*14.7+2.85005+0.08243*1.37775+0.23427*13.23137-0.08749*11.4
      , sigma = 0.58493-0.08516*14.7-0.17233-0.01741*11.4+0.10119*1.37775+0.08296*13.2-0.11659*13.23137
      , nu = 2.4123-0.4389*14.7-2.0214-3.7861*1.37775-1.8932*11.4+0.671*13.2+0.4136*13.23137, lower.tail = T, 
      log.p = F)






#for checking 
ECtest<-read.delim("E:/ECMWFdata/ECME_2019/Oct/ECME_NNL_2019100112/ECME_NNL_201910011200_NL096_LC",sep = "", header = F) #读取xlsx文件


#####remove the boundary############################################################################
rvi_table = read.delim("E:/KNMI_DATA/matchdata.txt", sep = "\t")
ra = raster("D:/ITC/KNMI/radar_data/09-20211013T195924Z-001/09/RAD_NL25_RAC_MFBS_01H_202009010200.nc")
ra = as.data.frame(ra)
ra$raid = c(1:535500)
ra_nl = na.omit(ra)
nl.merge = merge(rvi_table,ra_nl,by.x = "rID", by.y = "raid",all = TRUE)
nl.merge1 = merge(rvi_table,ra_nl,by.x = "rID", by.y = "raid")

t1 =  as.data.frame(table(nl.merge$df2id))
t2 = as.data.frame(table(nl.merge1$df2id))
colnames(t1) = c("dfID","all")
colnames(t2) = c("dfID","non-NA")
t = merge(t1,t2,by.x = "dfID",by.y = "dfID", all = TRUE)
t[is.na(t)] <- 0
t$percentage = t$`non-NA`/t$all
table(t$percentage == 1)
table(t$percentage >0.8)

table(t$`non-NA` >300)
table(t$`non-NA` >250)




















#####model the data###############################################


test = gamlss(formula = R4final1~EC4meanfinal1, sigma.formula = log(EC4sdfinal1),fimily=BCT,data=final,method=mixed(1,20))

test = ZAGA(mu.link = "log", sigma.link = "log", nu.link = "logit"))







objective_fun <- function(par, ensmean, ensvar, obs){
  m <- cbind(1, ensmean) %*% par[1:2]
  ssq_tmp <- cbind(1, ensvar) %*% par[3:4]
  if(any(ssq_tmp < 0)){
    return(999999)
  } else{
    s <- sqrt(ssq_tmp)
    return(sum(crps_norm(y = obs, location = m, scale = s)))
  }
}

gradfun_wrapper <- function(par, obs, ensmean, ensvar){
  loc <- cbind(1, ensmean) %*% par[1:2]
  sc <- sqrt(cbind(1, ensvar) %*% par[3:4])
  dcrps_dtheta <- gradcrps_norm(y = obs, location = loc, scale = sc) 
  out1 <- dcrps_dtheta[,1] %*% cbind(1, ensmean)
  out2 <- dcrps_dtheta[,2] %*% 
    cbind(1/(2*sqrt(par[3]+par[4]*ensvar)), 
          ensvar/(2*sqrt(par[3]+par[4]*ensvar)))
  return(as.numeric(cbind(out1,out2)))
}



postproc_global <- function(vdate, train_length){
  
  # determine training set
  train_end <- vdate - days(2)
  train_start <- train_end - days(train_length - 1)
  
  data_train <- subset(data, date >= train_start & date <= train_end)
  
  # remove incomplete cases (= NA obs or fc)
  data_train <- data_train[complete.cases(data_train), ]
  
  # determine optimal EMOS coefficients a,b,c,d using minimum CRPS estimation
  optim_out <- optim(par = c(1,1,1,1), 
                     fn = objective_fun,
                     gr = gradfun_wrapper,
                     ensmean = final$EC4meanfinal1, 
                     ensvar = (final$EC4sdfinal1)^2, 
                     obs = final$R4final1,
                     method = "BFGS")
  
  # check convergence of the numerical optimization
  if(optim_out$convergence != 0){
    message("numerical optimization did not converge")
  }
  
  # return optimal parameters
  return(optim_out$par)
}

setwd = "E:/RadarData/2018"

first_category_name = list.files("E01")
dir = paste("./11/",first_category_name,sep="")
n = length(dir) 


setwd('E:/ECMWFdata/ECME_2021')
# data = read.delim("E:/ECMWFdata/ECME_2019/Dec/ECME_NNL_2019121512/ECME_NNL_201912151200_NL084_LC",sep = "", header = F)
first_category_name = list.files("Jan")


n = length()
dates <- as.POSIXct(ncvar_get(nc, "time"), origin = "1970-01-01 00:00", tz = "UTC")
#######create a list for date#############################################################################################3

data1= as.data.frame(matrix(nrow = 4343, ncol = 1))
data1[,1] = as.POSIXct(data1[,1])

for (i in 1:4343) {
dates<- as.POSIXct((i-2)*3600, origin = "2018-11-01 00:00", tz = "UTC")
data1[i,] = dates
}

data2= as.data.frame(matrix(nrow = 5111, ncol = 1))
data2[,1] = as.POSIXct(data2[,1])

for (i in 1:5111) {
  dates<- as.POSIXct((i-2)*3600, origin = "2019-10-01 00:00", tz = "UTC")
  data2[i,] = dates
}

data3= as.data.frame(matrix(nrow = 5087, ncol = 1))
data3[,1] = as.POSIXct(data3[,1])

for (i in 1:5087) {
  dates<- as.POSIXct((i-2)*3600, origin = "2020-10-01 00:00", tz = "UTC")
  data3[i,] = dates
}

dates = rbind.data.frame(data1,data2,data3)
dates = dates[]

# example for one day
# tt <- as.Date("2016-01-01 00:00 UTC")
# m <- 50
# par_out <- postproc_global(tt, m)
# par_out
# data_eval <- subset(data, date == tt)
# loc <- c(cbind(1, data_eval$t2m_mean) %*% par_out[1:2])
# scsquared_tmp <- c(cbind(1, data_eval$t2m_var) %*% par_out[3:4])
# if(any(scsquared_tmp <= 0)){
#   print("negative scale, taking absolute value")
#   sc <- sqrt(abs(scsquared_tmp))
# } else{
#   sc <- sqrt(scsquared_tmp)
# }
# summary(crps_norm(y = data_eval$obs, mean = loc, sd = sc))





postproc_global <- function(vdate, train_length){
  
  # determine training set
  train_end <- vdate - days(2)
  train_start <- train_end - days(train_length - 1)
  
  data_train <- subset(data, date >= train_start & date <= train_end)
  
  # remove incomplete cases (= NA obs or fc)
  data_train <- data_train[complete.cases(data_train), ]
  
  # determine optimal EMOS coefficients a,b,c,d using minimum CRPS estimation
  optim_out <- optim(par = c(1,1,1,1), 
                     fn = objective_fun,
                     gr = gradfun_wrapper,
                     ensmean = final$EC4meanfinal1, 
                     ensvar = final$EC4sdfinal1, 
                     obs = final$R4final1,
                     method = "BFGS")
  
  
  
  
  
  
  # check convergence of the numerical optimization
  if(optim_out$convergence != 0){
    message("numerical optimization did not converge")
  }
  
  # return optimal parameters
  return(optim_out$par)
}



















#Reorganize data
ECmean <- ECmean %>%
  select(second_category, day, time, Rtime, grid, everything())
ECsd <- ECsd %>%
  select(second_category, day, time, Rtime, grid, everything())





R6h_1 = Rdata[,1:3]
EC6h_1 = cbind.data.frame(ECmean[1:242,1:3],ECsd[1:242,3])
colnames(EC6h_1) = c("first_category", "second_category","mean","sd")



dataa <- data[which(data$Dept_Code == "PEK" & substr(data$Plan_Dept_Time,7,7)==6 & substr(data$Plan_Dept_Time,9,9)==1
                    & substr(data$Plan_Dept_Time,10,10)==0),]
EC6h = ECmean[]




upscal_match = cbind.data.frame(radar_output[,1],radar_output[,5])
colnames(upscal_match) = c("rID","dfID")
write.table(upscal_match, "upscal_match.txt", sep = "\t", col.names = T, row.names = F)


#######################################################################################################        
  WDIR_ENS = WDIR_6hstep[3:53,]
    WDIR_ENS = as.data.frame(lapply(WDIR_ENS,as.numeric))
    for(k in 1:ncol(WDIR_ENS)){ 
      avg_var <- mean(WDIR_ENS[, k])
      avg_var_sin <- sin(avg_var*pi/180)
      avg_var_cos <- cos(avg_var*pi/180)
      median_var <- median(WDIR_ENS[,k])
      median_var_sin <- sin(median_var*pi/180)
      median_var_cos <- cos(median_var*pi/180)
      q10_var = quantile(WDIR_ENS[,K],probs=0.1)
      q10_var_sin <- sin(q10_var*pi/180)
      q10_var_cos <- cos(q10_var*pi/180)
      q90_var = quantile(WDIR_ENS[,K],probs=0.9)
      q90_var_sin <- sin(q90_var*pi/180)
      q90_var_cos <- cos(q90_var*pi/180)
#      new_data[i*j*9-8,]=WDIR_HRES
#      new_data[i*j*9-7,]=avg_var_sin
#      new_data[i*j*9-6,]=avg_var_cos
#      new_data[i*j*9-5,]=median_var_sin
#      new_data[i*j*9-4,]=median_var_cos
#      new_data[i*j*9-3,]=q10_var_sin
#      new_data[i*j*9-2,]=q10_var_cos
#      new_data[i*j*9-1,]=q90_var_sin
#      new_data[i*j*9,]=q90_var_cos
    }
  }
}
      
      
      new_1 = new_1[480:531,]
      new_1 = as.data.frame(lapply(new_1,as.numeric))
      new_1 = colMeans(new_1)
      names(new_1)<-c(1:84)
      new_1$third_category = 
        new_1$second_category<-substr(b[j],14,27)        #二级目录的名称是xlsx的文件名。
      new_1$first_category<-first_category_name[i]   #一级目录的名称是“文件夹名”
      rr_1<-rbind(merge_1,new_1)
    }
  }       

      
#added variables
13241
14034

20010
20257
00182   
      
      
    #WSPD (m/s)
    data_WSPD = data[161:213,]
    data_WSPD[,48:83]=trunc(data_WSPD[,48:83]/1000000)
    data_WSPD[2,64]=data_WSPD[2,64]/10^234
    data_WSPD = as.data.frame(lapply(data_WSPD,as.numeric))
    data_WSPD[2:53,]=data_WSPD[2:53,]/10
    WSPD_6hstep = data_WSPD[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
    #T2m (degrees Celcius)
    data_T2m = data[267:319,]
    data_T2m[,48:83]=trunc(data_T2m[,48:83]/1000000)
    data_T2m[2,64]=data_T2m[2,64]/10^234
    data_T2m = as.data.frame(lapply(data_T2m,as.numeric))
    data_T2m[2:53,]=data_T2m[2:53,]/10
    T2m_6hstep = data_T2m[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
    #DPT (degrees Celcius)
    data_DPT = data[320:372,]
    data_DPT[,48:83]=trunc(data_DPT[,48:83]/1000000)
    data_DPT[2,64]=data_DPT[2,64]/10^234
    data_DPT = as.data.frame(lapply(data_DPT,as.numeric))
    data_DPT[2:53,]=data_DPT[2:53,]/10
    DPT_6hstep = data_DPT[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
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
    #SNOWcum (mm)
    data_SNOWcum = data[532:584,]
    data_SNOWcum[,48:83]=trunc(data_SNOWcum[,48:83]/1000000)
    data_SNOWcum[2,64]=data_SNOWcum[2,64]/10^234
    data_SNOWcum = as.data.frame(lapply(data_SNOWcum,as.numeric))
    data_SNOWcum[2:53,]=data_SNOWcum[2:53,]/10
    SNOWcum_6hstep = data_SNOWcum[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
    #SNOW (mm)
    SNOW_6hstep = SNOWcum_6hstep[,1:(ncol(SNOWcum_6hstep)-1)]
    SNOW_6hstep = cbind(0,SNOW_6hstep)
    SNOWcum_6hstep = as.data.frame(lapply(SNOWcum_6hstep,as.numeric))
    SNOW_6hstep = as.data.frame(lapply(SNOW_6hstep,as.numeric))
    SNOW_6hstep = SNOWcum_6hstep-SNOW_6hstep
    #CC (%)
    data_CC = data[638:690,]
    data_CC[,48:83]=trunc(data_CC[,48:83]/1000000)
    data_CC[2,64]=data_CC[2,64]/10^234
    data_CC = as.data.frame(lapply(data_CC,as.numeric))
    CC_6hstep = data_CC[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
    #conv precip_cum 
    data_conRRcum = data[691:743,]
    data_conRRcum[,48:83]=trunc(data_conRRcum[,48:83]/1000000)
    data_conRRcum[2,64]=data_conRRcum[2,64]/10^234
    data_conRRcum = as.data.frame(lapply(data_conRRcum,as.numeric))
    data_conRRcum[2:53,]=data_conRRcum[2:53,]/10
    conRRcum_6hstep = data_conRRcum[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
    #conv precip (mm)
    conRR_6hstep = conRRcum_6hstep[,1:(ncol(conRRcum_6hstep)-1)]
    conRR_6hstep = cbind(0,conRR_6hstep)
    conRRcum_6hstep = as.data.frame(lapply(conRRcum_6hstep,as.numeric))
    conRR_6hstep = as.data.frame(lapply(conRR_6hstep,as.numeric))
    conRR_6hstep = conRRcum_6hstep-conRR_6hstep
    #CAPE
    data_cape = data[585:637,]
    data_cape[,48:83]=trunc(data_cape[,48:83]/1000000)
    data_cape[2,64]=data_cape[2,64]/10^234
    data_cape = as.data.frame(lapply(data_cape,as.numeric))
    cape_6hstep = data_cape[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
    #sunshine_duration (the data has some problems)
    #CAPE shear
    data_capeS = data[849:901,]  #For Y2018 and Y2019 Jan, Feb, Mar, Apr, Oct(1-5)
    
    data_capeS = data[797:849,]  #Start from Oct 6th
    
    data_capeS[,48:83]=trunc(data_capeS[,48:83]/1000000)
    data_capeS[2,64]=data_capeS[2,64]/10^234
    data_capeS = as.data.frame(lapply(data_capeS,as.numeric))
    capeS_6hstep = data_capeS[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
    #evaporation_ cum
    data_evacum = data[955:1007,] #For Y2018 and Y2019 Jan, Feb, Mar, Apr, Oct(1-5)
    data_evacum = data[903:955,] #Start from Oct 6th
    
    data_evacum[,48:83]=trunc(data_evacum[,48:83]/1000000)
    data_evacum[2,64]=data_evacum[2,64]/10^234
    data_evacum = as.data.frame(lapply(data_evacum,as.numeric))
    evacum_6hstep = data_evacum[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
    #evaporation (mm)
    eva_6hstep = evacum_6hstep[,1:(ncol(evacum_6hstep)-1)]
    eva_6hstep = cbind(0,eva_6hstep)
    evacum_6hstep = as.data.frame(lapply(evacum_6hstep,as.numeric))
    eva_6hstep = as.data.frame(lapply(eva_6hstep,as.numeric))
    eva_6hstep = evacum_6hstep-eva_6hstep   
    
    

        df = WDIR_6hstep[3:53,]
    p1 <- matrix(0, 2, ncol(df))
#    for (i in 1:nrow(df)) {
      for (j in 1:ncol(df)) {
        avg_var <- mean(df[, j])
        sd_var  <- sd(df[, j])
        p1[1,j] = avg_var
        p1[2,j] = sd_var

      }

    
    
    
    

    

        
#    WDIR_6hstep_sd=apply(WDIR_6hstep[3:53,], 2, sd)   
    

    
    
    
    
    



#  write.xlsx(merge_1,paste(dir[i],'/merge.xlsx',sep=''),row.names = F,col.names= F)






#For the data from the year 2020
setwd("E:/ECMWFdata/ECME_2020/ECME_NNL_2020110100")
data = read.delim("ECME_NNL_202011010000_NL012_LC",sep = "", header = F)






sin(90*pi/180)
cos(90*pi/180)


#    WDIR_6hstep_HRES = WDIR_6hstep[2,]
#WDIR_6hstep_mean = colMeans(WDIR_6hstep[3:53,])
#WDIR_6hstep_median=apply(WDIR_6hstep[3:53,], 2, median)



setwd("E:/ECMWFdata/ECME_NNL_2020110100")
all_file = list.files()
n = length(all_file)
member = 52
lead_time = 40
file = paste("E:/ECMWFdata/ECME_NNL_2021110100", all_file, sep = "/")

df = as.data.frame(matrix(nrow = n*40, ncol = 54))
df2 = as.data.frame(matrix(nrow =242, ncol = 3))




colnames(df)[1] = 'lead-time'
df$`lead-time` = rep(1:40, each = 242)

colnames(df)[2] = 'grid'
df$grid = rep(c(1:242), times = 40)

for (j in 1:n){
  
  data = read.delim(file[j],sep = "", header = F)
#data = read.delim("ECME_NNL_201912010000_NL012_LC",sep = "", header = F)
  data_rain = data[c(480:531),c(1:40)]
  
  data_rain_1 = data_rain[,1:(ncol(data_rain)-1)]
  data_rain_1 = cbind(0,data_rain_1)
  
  data_rain = as.data.frame(lapply(data_rain,as.numeric))
  data_rain_1 = as.data.frame(lapply(data_rain_1,as.numeric))
  
  data_rain = data_rain-data_rain_1
  for (i in 1:40){
    df[(i-1)*242+j,c(3:54)] = data_rain[,i] }
  
}
for (k in 1:n){
  data = read.delim(file[k],sep = "", header = F)
  
  df2[k,] = data[1,1:3]
  
}






library(e1071)
library(gamlss.dist)
library(ddplyr)
mean(data)
sd(data)
skewness(data)
kurtosis(data)
quantile(data,probs = seq(0,1,0.25))





test <- test %>% dplyr::mutate(
  q025 = gamlss.dist::qZAGA(rep(0.025, n()), mu = mu, sigma = 1 / sqrt(sigma), nu = pzero),
  q975 = gamlss.dist::qZAGA(rep(0.975, n()), mu = mu, sigma = 1 / sqrt(sigma), nu = pzero),
  mean = mu)

#different lead time should have different model; but different grid can have the same model (spatially)
#high_resolution can be a separate variable; and check other variables at the same time
