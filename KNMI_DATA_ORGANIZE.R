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




saveRDS(object = df, file = "ECME_NNL_2020110100.rds")
tf2 = readRDS(file = "ECME_NNL_2020110100.rds")

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
ra = raster("D:/ITC/KNMI/radar_data/09-20211013T195924Z-001/09/RAD_NL25_RAC_MFBS_01H_202009010200.nc")
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
new_ra = aggregate(ra, fact=3, expand=FALSE, fun=mean, na.rm=TRUE)
new_ra1 = aggregate(new_ra, fact=3, expand=FALSE, fun=max, na.rm=TRUE)

lat = init(ra,'y')
lat = as.data.frame(lat)


#max RA 9*9
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
write.table(upscal_match_9km, "upscal_match_9km.txt", sep = "\t", col.names = T, row.names = F)

ggplot(radar_output) +
  geom_point(aes(x = rLon, y= rLat, color = image1_image_data)) 





setwd("E:/RadarData/2021")
first_category_name = list.files("03")
dir = paste("E:/RadarData/2021/03/",first_category_name,sep="")
n = length(dir) 

#radar_combine = as.data.frame(matrix(0,nrow = 599, ncol = n))
date = as.data.frame(matrix(0,nrow = 1, ncol = n))
#names(radar_combine) = c(1:n)

for(i in 1:n){
#  ra<-raster(dir[i]) 
  rb<-nc_open(dir[i])
#  ra[397,369] = ra[397,370]= ra[445,352] = ra[445,353] = ra[435,309] = ra[435,310] =
#    ra[436,309] = ra[436,310] = ra[504,281] = ra[504,282] = ra[504,283] = ra[505,282] =
#    ra[505,283] = ra[449,289] = ra[449,290] = ra[449,291] = ra[449,292] = ra[450,288] = 
#    ra[450,289] = ra[450,290] = ra[450,291] = ra[450,293] = ra[450,294] = ra[451,288] = 
#    ra[451,289] = ra[451,294] = ra[452,289] = ra[452,290] = ra[452,291] = ra[452,292] = 
#    ra[453,289] = ra[453,290] = ra[453,291] = ra[453,292] = ra[326,446] = ra[326,447] = 
#    ra[327,446] = ra[327,447] = ra[327,448] = ra[328,446] = ra[328,447] = ra[328,448] = 
#    ra[330,447] = NA
  
  
#  new_ra = aggregate(ra, fact=3, expand=FALSE, fun=mean, na.rm=TRUE)
#  new_ra1 = aggregate(new_ra, fact=3, expand=FALSE, fun=max, na.rm=TRUE)
  
#  data_matrix <- rasterToPoints(new_ra1)
#  data_matrix= as.data.frame(data_matrix)
#  output = data_matrix$image1_image_data
#  radar_combine[,i]=output
  dates <- as.POSIXct(ncvar_get(rb, "time"), origin = "2000-01-01 01:00", tz = "UTC")
  date[,i]=dates
}



new_df = radar_combine

f <- matrix(0, dim(new_df)[1], dim(new_df)[2] %/% 6 )
for (j in seq(1, (dim(new_df)[2] - dim(new_df)[2] %% 6), by = 6)) {
  f[, (j - 1) / 6 + 1] <- rowSums(new_df[, j:(j + 5)])
}

f = as.data.frame(f)
colnames(f) = 1: dim(f)[2]
setwd("E:/RadarData/radar_output")

write.table(f, "2021_04_9km.txt", sep = "\t", col.names = T, row.names = F)





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



####(for_RR)###############################################################################################    
library(e1071)
setwd('E:/ECMWFdata/ECME_2021')
# data = read.delim("E:/ECMWFdata/ECME_2019/Dec/ECME_NNL_2019121512/ECME_NNL_201912151200_NL084_LC",sep = "", header = F)
first_category_name = list.files("Jan")
dir = paste("./Jan/",first_category_name,sep="")
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

#####(For Snow)#################################################################################################
library(e1071)
setwd('E:/ECMWFdata/ECME_2021')
# data = read.delim("E:/ECMWFdata/ECME_2019/Dec/ECME_NNL_2019121512/ECME_NNL_201912151200_NL084_LC",sep = "", header = F)
first_category_name = list.files("Apr")
dir = paste("./Apr/",first_category_name,sep="")
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
    data_RRcum = data[532:584,]
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
write.table(merge_HRES, "Snow_6h_HRES_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_avg, "Snow_6h_mean_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_sd, "Snow_6h_sd_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_median, "Snow_6h_median_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q10, "Snow_6h_q10_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_q90, "Snow_6h_q90_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_p0, "Snow_6h_p0_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_p1, "Snow_6h_p1_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_p3, "Snow_6h_p3_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_p5, "Snow_6h_p5_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_p10, "Snow_6h_p10_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_p20, "Snow_6h_p20_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_skewness, "Snow_6h_skewness_2018.txt", sep = "\t", col.names = T, row.names = F)
write.table(merge_kurtosis, "Snow_6h_kurtosis_2018.txt", sep = "\t", col.names = T, row.names = F)



######Temporally match datasets######################################################################
#Read forecast data and radar data
library(dplyr)
ECmean201811 = read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_mean_2018_Nov.txt",sep = "", header = T)
ECmean201811 = ECmean201811[-1,]
ECmean201812 = read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_mean_2018_Dec.txt",sep = "", header = T)
ECmean201812 = ECmean201812[-1,]
ECmean201901 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_mean_2019_Jan.txt",sep = "", header = T)
ECmean201901 = ECmean201901[-1,]
ECmean201902 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_mean_2019_Feb.txt",sep = "", header = T)
ECmean201902 = ECmean201902[-1,]
ECmean201903 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_mean_2019_Mar.txt",sep = "", header = T)
ECmean201903 = ECmean201903[-1,]
ECmean201904 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_mean_2019_Apr.txt",sep = "", header = T)
ECmean201904 = ECmean201904[-1,]
ECmean201910 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_mean_2019_Oct.txt",sep = "", header = T)
ECmean201910 = ECmean201910[-1,]
ECmean201911 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_mean_2019_Nov.txt",sep = "", header = T)
ECmean201911 = ECmean201911[-1,]
ECmean201912 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_mean_2019_Dec.txt",sep = "", header = T)
ECmean201912 = ECmean201912[-1,]
ECmean202001 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_mean_2020_Jan.txt",sep = "", header = T)
ECmean202001 = ECmean202001[-1,]
ECmean202002 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_mean_2020_Feb.txt",sep = "", header = T)
ECmean202002 = ECmean202002[-1,]
ECmean202003 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_mean_2020_Mar.txt",sep = "", header = T)
ECmean202003 = ECmean202003[-1,]
ECmean202004 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_mean_2020_Apr.txt",sep = "", header = T)
ECmean202004 = ECmean202004[-1,]
ECmean202010 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_mean_2020_Oct.txt",sep = "", header = T)
ECmean202010 = ECmean202010[-1,]
ECmean202011 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_mean_2020_Nov.txt",sep = "", header = T)
ECmean202011 = ECmean202011[-1,]
ECmean202012 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_mean_2020_Dec.txt",sep = "", header = T)
ECmean202012 = ECmean202012[-1,]
ECmean202101 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_mean_2021_Jan.txt",sep = "", header = T)
ECmean202101 = ECmean202101[-1,]
ECmean202102 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_mean_2021_Feb.txt",sep = "", header = T)
ECmean202102 = ECmean202102[-1,]
ECmean202103 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_mean_2021_Mar.txt",sep = "", header = T)
ECmean202103 = ECmean202103[-1,]
ECmean202104 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_mean_2021_Apr.txt",sep = "", header = T)
ECmean202104 = ECmean202104[-1,]







ECsd201811 = read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_sd_2018_Nov.txt",sep = "", header = T)
ECsd201811 = ECsd201811[-1,]
ECsd201812 = read.table("E:/ECMWFdata/pre-process/Y2018/RR_6h_sd_2018_Dec.txt",sep = "", header = T)
ECsd201812 = ECsd201812[-1,]
ECsd201901 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_sd_2019_Jan.txt",sep = "", header = T)
ECsd201901 = ECsd201901[-1,]
ECsd201902 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_sd_2019_Feb.txt",sep = "", header = T)
ECsd201902 = ECsd201902[-1,]
ECsd201903 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_sd_2019_Mar.txt",sep = "", header = T)
ECsd201903 = ECsd201903[-1,]
ECsd201904 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_sd_2019_Apr.txt",sep = "", header = T)
ECsd201904 = ECsd201904[-1,]
ECsd201910 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_sd_2019_Oct.txt",sep = "", header = T)
ECsd201910 = ECsd201910[-1,]
ECsd201911 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_sd_2019_Nov.txt",sep = "", header = T)
ECsd201911 = ECsd201911[-1,]
ECsd201912 = read.table("E:/ECMWFdata/pre-process/Y2019/RR_6h_sd_2019_Dec.txt",sep = "", header = T)
ECsd201912 = ECsd201912[-1,]
ECsd202001 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_sd_2020_Jan.txt",sep = "", header = T)
ECsd202001 = ECsd202001[-1,]
ECsd202002 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_sd_2020_Feb.txt",sep = "", header = T)
ECsd202002 = ECsd202002[-1,]
ECsd202003 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_sd_2020_Mar.txt",sep = "", header = T)
ECsd202003 = ECsd202003[-1,]
ECsd202004 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_sd_2020_Apr.txt",sep = "", header = T)
ECsd202004 = ECsd202004[-1,]
ECsd202010 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_sd_2020_Oct.txt",sep = "", header = T)
ECsd202010 = ECsd202010[-1,]
ECsd202011 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_sd_2020_Nov.txt",sep = "", header = T)
ECsd202011 = ECsd202011[-1,]
ECsd202012 = read.table("E:/ECMWFdata/pre-process/Y2020/RR_6h_sd_2020_Dec.txt",sep = "", header = T)
ECsd202012 = ECsd202012[-1,]
ECsd202101 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_sd_2021_Jan.txt",sep = "", header = T)
ECsd202101 = ECsd202101[-1,]
ECsd202102 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_sd_2021_Feb.txt",sep = "", header = T)
ECsd202102 = ECsd202102[-1,]
ECsd202103 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_sd_2021_Mar.txt",sep = "", header = T)
ECsd202103 = ECsd202103[-1,]
ECsd202104 = read.table("E:/ECMWFdata/pre-process/Y2021/RR_6h_sd_2021_Apr.txt",sep = "", header = T)
ECsd202104 = ECsd202104[-1,]






#Rdata201811 =  read.table("E:/RadarData/radar_output/2018_11.txt",sep = "", header = T)
Rdata = cbind.data.frame(read.table("E:/output/radar_output/2018_11_9km.txt",sep = "", header = T),
                             read.table("E:/output/radar_output/2018_12_9km.txt",sep = "", header = T),
                             read.table("E:/output/radar_output/2019_01_9km.txt",sep = "", header = T),
                             read.table("E:/output/radar_output/2019_02_9km.txt",sep = "", header = T),
                             read.table("E:/output/radar_output/2019_03_9km.txt",sep = "", header = T),
                             read.table("E:/output/radar_output/2019_04_9km.txt",sep = "", header = T),
                             read.table("E:/output/radar_output/2019_10_9km.txt",sep = "", header = T),
                             read.table("E:/output/radar_output/2019_11_9km.txt",sep = "", header = T),
                             read.table("E:/output/radar_output/2019_12_9km.txt",sep = "", header = T),
                             read.table("E:/output/radar_output/2020_01_9km.txt",sep = "", header = T),
                             read.table("E:/output/radar_output/2020_02_9km.txt",sep = "", header = T),
                             read.table("E:/output/radar_output/2020_03_9km.txt",sep = "", header = T),
                             read.table("E:/output/radar_output/2020_04_9km.txt",sep = "", header = T),
                             read.table("E:/output/radar_output/2020_10_9km.txt",sep = "", header = T),
                             read.table("E:/output/radar_output/2020_11_9km.txt",sep = "", header = T),
                             read.table("E:/output/radar_output/2020_12_9km.txt",sep = "", header = T),
                             read.table("E:/output/radar_output/2021_01_9km.txt",sep = "", header = T),
                             read.table("E:/output/radar_output/2021_02_9km.txt",sep = "", header = T),
                             read.table("E:/output/radar_output/2021_03_9km.txt",sep = "", header = T),
                             read.table("E:/output/radar_output/2021_04_9km.txt",sep = "", header = T))
                             
#new = read.table("E:/output/time_information.txt", sep = "\t", header = T)
#new = t(new)
#colnames(new) = c(1:2424)
#colnames(Rdata) = c(1:2424)
#Rdata = rbind.data.frame(new,Rdata)
                             
#write.table(Rdata, "time_information.txt", sep = "\t", col.names = T, row.names = F)

                             
                             
                             
library(lubridate)                             

upscal_match = read.table("E:/KNMI_DATA/upscal_match_9km.txt", sep = "", header = T)
colnames(Rdata) = c(1:dim(Rdata)[2])
#Match data temporally
Rdata$rID = upscal_match$rID
Rdata$dfID = upscal_match$dfID

#add time information for EC data
ECmean = rbind.data.frame(ECmean201811, ECmean201812, ECmean201901, ECmean201902, ECmean201903, ECmean201904, ECmean201910, ECmean201911,
                          ECmean201912, ECmean202001, ECmean202002, ECmean202003, ECmean202004, ECmean202010, ECmean202011, ECmean202012,
                          ECmean202101, ECmean202102, ECmean202103, ECmean202104)
ECsd = rbind.data.frame(ECsd201811, ECsd201812, ECsd201901, ECsd201902, ECsd201903, ECsd201904, ECsd201910, ECsd201911, ECsd201912,
                        ECsd202001, ECsd202002, ECsd202003, ECsd202004, ECsd202010, ECsd202011, ECsd202012, ECsd202101, ECsd202102,
                        ECsd202103, ECsd202104)

rm(ECmean201811,ECmean201812,ECmean201901,ECmean201902,ECmean201903,ECmean201904,ECmean201910,ECmean201911,
   ECmean201912,ECmean202001,ECmean202002,ECmean202003,ECmean202004,ECmean202010,ECmean202011,ECmean202012,
   ECmean202101,ECmean202102,ECmean202103,ECmean202104,ECsd201811,ECsd201812,ECsd201901,ECsd201902,ECsd201903,
   ECsd201904,ECsd201910,ECsd201911,ECsd201912,ECsd202001,ECsd202002,ECsd202003,ECsd202004,ECsd202010,
   ECsd202011,ECsd202012,ECsd202101,ECsd202102,ECsd202103,ECsd202104)    


ECmean$grid = as.numeric(substr(ECmean$second_category,12,15))
ECsd$grid = as.numeric(substr(ECsd$second_category,12,15))


ECmean$Rtime = as.POSIXct("2000-1-1 01:00:00", tz = "UTC")
year(ECmean$Rtime) = as.numeric(substr(ECmean$first_category,10,13))
month(ECmean$Rtime) = as.numeric(substr(ECmean$first_category,14,15))
mday(ECmean$Rtime) = as.numeric(substr(ECmean$first_category,16,17))
hour(ECmean$Rtime) = as.numeric(substr(ECmean$first_category,18,19))
minute(ECmean$Rtime) = 00
second(ECmean$Rtime) = 00

ECsd$Rtime = as.POSIXct("2000-1-1 01:00:00", tz = "UTC")
year(ECsd$Rtime) = as.numeric(substr(ECsd$first_category,10,13))
month(ECsd$Rtime) = as.numeric(substr(ECsd$first_category,14,15))
mday(ECsd$Rtime) = as.numeric(substr(ECsd$first_category,16,17))
hour(ECsd$Rtime) = as.numeric(substr(ECsd$first_category,18,19))
minute(ECsd$Rtime) = 00
second(ECsd$Rtime) = 00









   #outputD = Rdata[,c(4,2425,2426)]

   

   inputmean = as.data.frame(ECmean[,c(4,61:65)])
   inputsd = as.data.frame(ECsd[,c(4,61:65)])



   Rdata1 = Rdata[,-c(2425,2426)]
#   Rdata1 = Rdata1[,-c(1:3)]
   
   EC4meanfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 599, ncol = 1))
   EC4sdfinal <- as.data.frame(matrix(upscal_match$rID,nrow = 599, ncol = 1))
   R4final <- as.data.frame(matrix(upscal_match$rID,nrow = 599, ncol = 1))
   
#make df for training   
   for (i in unique(inputmean$Rtime)){

     subdata = cbind.data.frame((Rdata1[,i+4]),Rdata[,2425],Rdata[,2426])
     colnames(subdata) = c("Rdata","RID","dfID")
     submean = inputmean[which(inputmean$Rtime == i),]
     subsd = inputsd[which(inputsd$Rtime == i),]
     subpre = cbind.data.frame(submean[,1],subsd)
     colnames(subpre)[1] = "submean"
     colnames(subpre)[2] = "subsd"
     subpre = subpre[,c(1,2,6,7)]
     match = left_join(subdata,subpre,by = c("dfID"="grid"))
     EC4Rmean = match[,4]
     EC4sd = match[,5]
     R4match = match[,1]
     EC4meanfinal=cbind.data.frame(EC4meanfinal,EC4Rmean)
     EC4sdfinal = cbind.data.frame(EC4sdfinal,EC4sd)
     R4final = cbind.data.frame(R4final,R4match)

   }
list =   unique(inputmean$Rtime) 

rm(ECmean201811,ECmean201812,ECmean201901,ECmean201902,ECmean201903,ECmean201904,ECmean201910,ECmean201911,
   ECmean201912,ECmean202001,ECmean202002,ECmean202003,ECmean202004,ECmean202010,ECmean202011,ECmean202012,
   ECmean202101,ECmean202102,ECmean202103,ECmean202104,ECsd201811,ECsd201812,ECsd201901,ECsd201902,ECsd201903,
   ECsd201904,ECsd201910,ECsd201911,ECsd201912,ECsd202001,ECsd202002,ECsd202003,ECsd202004,ECsd202010,
   ECsd202011,ECsd202012,ECsd202101,ECsd202102,ECsd202103,ECsd202104)     
rm(ECmean,ECsd,inputmean,inputsd,Rdata,Rdata1,subdata,submean,subpre,subsd,upscal_match)
   

EC4meanfinal1 = EC4meanfinal[,-1]
EC4meanfinal1 = as.matrix(EC4meanfinal1)
dim(EC4meanfinal1) <- c(599*1120,1)
EC4sdfinal1 = EC4sdfinal[,-1]
EC4sdfinal1 = as.matrix(EC4sdfinal1)
dim(EC4sdfinal1) <- c(599*1120,1)
R4final1 = R4final[,-1]
R4final1 = as.matrix(R4final1)
dim(R4final1) <- c(599*1120,1)

41*1120=45920
74*1120=82880
77*1120=86240
89*1120=99680
250/396=0.631

final = cbind.data.frame(EC4meanfinal1,EC4sdfinal1,R4final1)
final[which(final$EC4meanfinal1 <0 ),1]<- 0

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




















#####model with mean and variance###############################################
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
    CC_6hstep = data_CC[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,49:84)]
 
    avg_WDIR = mean()
    
    
    
    
    

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
