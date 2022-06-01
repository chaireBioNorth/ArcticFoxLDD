
###############################################
#####         Create LDD tables          #####
###############################################

# by Richard Gravel


###
library(readxl)
library(proj4)
library(dplyr)
library(rgdal)
library(raster)
library(sp)
library(patchwork)
library(hrbrthemes)
library(zoo)
library(ggplot2)
library(spatstat)
library(maptools)
library(rgeos)
library(tidyr)
library(geosphere)

#Dataset
#Download database (transience only AND all locs)
setwd( "")
fox_transience <- read.delim("argos_transience.txt", header = T)
fox <- read.delim("argos_filtrees_byday_2007_2021.txt", header = T)

#Create objects for projections
wgs1984.proj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
UTMZ17 <- CRS("+proj=utm +zone=17 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
polar <- CRS("+proj=stere +lat_0=90 +lat_ts=71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

#TABLE: cols
setwd("/Users/richardgravel/Documents/Manuscrit/data")
tab_Collar <- read.table("tab_LDD.txt", header = T)
cols <- tab_Collar
colnames(cols) <- c("Fox_ID", "Platform_ID","start_tracking","start_disp","end_disp","end_tracking","end_status","life_stage", "disp_status","sex","migration_merged")

#List of dispersal foxes
id_disp <-   levels(as.factor(as.character(cols[,"Platform_ID"])))


#Subset data base (keeping dispersal foxes)
fox <- fox[fox$Platform_ID %in% id_disp,]

fox <- mutate_at(fox, vars(Fox_ID, Platform_ID, Loc_quality), list(factor))
str(fox)

fox$Loc_date <- as.POSIXct(strptime(fox$Loc_date, format = "%Y-%m-%d %H:%M:%S", tz = "GMT"))

table(droplevels(fox)$Fox_ID)
table(droplevels(fox)$Platform_ID)


#Project the data
fox_sp <- SpatialPointsDataFrame(fox[, 3:4], data = fox[, -c(3:4)], proj4string = wgs1984.proj)
fox_sp <- spTransform(fox_sp, polar)



setwd("/Digital_Chart_of_the_World")

#Create object for database
xy <- fox_sp

#Import Bylot outline and change CRS
World <- rgdal::readOGR("Digital_Chart_of_the_World.shp")
World_Polar <- spTransform(World, polar)

#Polygon into SpatialPolygons format
world_pol <- World_Polar[,1]@polygons
world <- SpatialPolygons(world_pol, proj4string = polar)

# Clip Bylot on points for all foxes
clipped <- over(xy, world)

# Points IN Bylot
xy_in <- xy[which(!is.na(clipped)), ]

# Points OUT Bylot
xy_out <- xy[which(is.na(clipped)), ]

dat <- as.data.frame(matrix(nrow = length(xy_out), ncol = 1))
In0_Out1World <- rep("1", nrow(xy_out))
comp <- as.data.frame(xy_out)
yx2 <- cbind(comp, In0_Out1World)

dat <- as.data.frame(matrix(nrow = length(xy_in), ncol = 1))
In0_Out1World <- rep("0", nrow(xy_in))
comp <- as.data.frame(xy_in)
xy_in2 <- cbind(comp, In0_Out1World)

tab1 <- rbind(yx2, xy_in2)
tab1 <- tab1[order(tab1[, "Fox_ID"],tab1[, "Loc_date"]), ]


#Verification of the classification In and Out (to do in UTM so i can open it on QGIS)
#setwd("/Users/richardgravel/Desktop/Fox/masterfox/GIS")
#INworld <- subset(tab1, In0_Out1World==0)
#INworld <- SpatialPointsDataFrame(INworld[, 5:6], data = INworld[, -c(5:6)], proj4string = polar)
#Inworld <- spTransform(INworld, UTMZ17)
#shapefile(INworld, filename = "INworld", overwrite= TRUE)

#OUTworld <- subset(tab1, In0_Out1World==1)
#OUTworld <- SpatialPointsDataFrame(OUTworld[, 5:6], data = OUTworld[, -c(5:6)], proj4string = polar)
#OUTworld <- spTransform(OUTworld, UTMZ17)
#shapefile(OUTworld, filename = "OUTworld", overwrite= TRUE)


#TABLE: tab (distances)

#join tab1$In0_Out1World to fox_sp

join.fox <- left_join(fox, tab1, by=c("Fox_ID","Loc_date"))
fox <- as.data.frame(c(join.fox[1:6], join.fox[11]))
colnames(fox) <- c("Fox_ID","Platform_ID","Longitude","Latitude","Loc_date","Loc_quality","In0_Out1_Land")

fox_sp <- SpatialPointsDataFrame(fox[, 3:4], data = fox[, -c(3:4)], proj4string = wgs1984.proj)
fox_sp <- spTransform(fox_sp, polar)


#Create new empty table 
tab2 <- data.frame()

#Loop by ID

for (i in 1 : length(id_disp)){
  
  XY <- fox_sp[fox_sp@data$Platform_ID == id_disp[[i]],]
  
  
  #Mean of the 14 firsts locs as center
  sub <- (fox_sp[fox_sp$Platform_ID%in%id_disp[[i]], ][2:15, ])
  cen <- c(mean(coordinates(sub)[,1]), mean(coordinates(sub)[,2]))
  cen_sp_UTM <- SpatialPoints(data.frame(cen[1], cen[2]), proj4string =polar)
  sub_first <- sub[1, ]
  
  center <- as.data.frame(matrix(nrow = 1, ncol = 7))
  colnames(center) <- c("Fox_ID","Platform_ID","Longitude","Latitude","Loc_date","Loc_quality","In0_Out1_Land")
  center[,1:2] <- sub_first@data[1:2]
  center[,5:7] <- NA
  center[,3:4] <- cen
  
  center <- SpatialPointsDataFrame(coords = center[,3:4], data = center[,c(1:2,5:7)], proj4string = polar)
  
  # Loc's distance form center 
  #===================================
  
  # Calculation of minimal Eucledean distances from center
  
  distance <- rep(NA,nrow(XY@coords))
  ppp_points <- as(as(XY, "SpatialPoints"), "ppp")  # convert point format
  psp_center <- as.ppp.SpatialPointsDataFrame(center) # convert center in ppp format
  
  for (j in 1 : length(distance)){
    
    #Distance of loc j from centre
    Dist <- nncross(ppp_points[j,], psp_center)$dist    
    distance[j] <- Dist        
  }
  
  #Create a matrix of distances for one individual
  #Matrix with distances of points from center
  mat <- as.data.frame(matrix(nrow = length(XY), ncol = 8))
  colnames(mat) <- c("Fox_ID","Platform_ID","Longitude","Latitude","Loc_date","Loc_quality","In0_Out1_Land","Distance_Km")
  mat[,1:2] <- XY@data[1:2]
  mat[,5:7] <- XY@data[3:5]
  mat[,3:4] <- XY@coords
  mat[,8] <- distance/1000
  
  res <- mat
  
  
  #Order matrix by date
  ord <- order(res[,"Loc_date"])
  res <- res[ord,]
  
  #Stock into dataframe
  tab2 <- rbind(tab2, res)
}

tab2 <- tab2[order(tab2[,"Fox_ID"],tab2[,"Loc_date"]),]




join.fox2 <- inner_join(fox, tab2, by=c("Fox_ID","Loc_date"))
fox2 <- as.data.frame(c(join.fox2[1:7], join.fox2[13]))
colnames(fox2) <- c("Fox_ID","Platform_ID","Longitude","Latitude","Loc_date","Loc_quality","In0_Out1_Land","Distance_Km")

fox2_sp <- SpatialPointsDataFrame(fox2[, 3:4], data = fox2[, -c(3:4)], proj4string = wgs1984.proj)
fox2_sp <- spTransform(fox2_sp, polar)



#Daily table of dispersal foxes (Daily rate, NSD/MSD, In0_Out1_Land)

#Create empty dataframe
table_daily_disp <- data.frame()

for (i in 1 : length(id_disp)){
  
  #subsets for tables
  xx <- fox2_sp[fox2_sp$Platform_ID == id_disp[[i]],]
  yy <- cols[cols$Platform_ID== id_disp[[i]],]
  
  
  #fox_ID, sex and life_stage
  fox_ID <- as.character(xx$Fox_ID[1])
  sex <- as.character(yy$sex)
  life_stage <- as.character(yy$life_stage)

  
  #FOR DAILY RATE 
  
  #Create a time lag matrix
  time <- as.POSIXct(xx$Loc_date,"%Y-%m-%d",tz = "GMT")
  time <- substr(as.character(xx$Loc_date),1,10)
  time_lag <- as.numeric(NULL)
  
  for (j in 1 : length(time)){
    Lag <-  abs(difftime(time[j], time[j-1], units='days'))
    time_lag[j] <- Lag  
  }
  
  
  #Create a distance lag matrix (km)
  distance <- rep(NA,nrow(xx@coords))
  ppp_points <- as(as(xx, "SpatialPoints"), "ppp")  # convert point format
  
  for (j in 1 : length(distance)){
    Dist <- nncross(ppp_points[j,], ppp_points[j-1,])$dist    
    distance[j] <- Dist/1000     
  }
  distance[1] <-  0
  
  
  #Create a daily rate matrix + other importants columns
  rate <- as.data.frame(matrix(nrow = length(distance), ncol = 6))
  colnames(rate) <- c("Fox_ID","tracking_day","Loc_date","daily_rate","In0_Out1_Land","Distance_Km")
  rate[,1] <- fox_ID
  rate[,2] <- NA
  rate[,3] <- as.POSIXct(strptime(xx$Loc_date,"%Y-%m-%d",tz = "GMT"))
  rate[,4] <- distance/time_lag
  rate[,5] <- xx$In0_Out1_Land
  rate[,6] <- xx$Distance_Km

  
  
  #Create dataframe with all date (including when location is missing)
  xx$Loc_date <- as.POSIXct(strptime(xx$Loc_date,"%Y-%m-%d",tz = "GMT"))
  ts <- seq.POSIXt(xx$Loc_date[1], xx$Loc_date[nrow(xx)], by="day")
  df <- data.frame(Loc_date=ts)
  
  
  #full.rate <- full_join(rate, full_date)
  full_rate <- full_join(rate, df)
  full_rate <- full_rate[order(full_rate[,"Loc_date"]),]
  full_rate[,1] <- fox_ID
  full_rate[,2] <- seq.int(from=0,(nrow(full_rate)-1))
  
  #Fill missing rate value (giving the rate value that succeded the missing loc)
  full_rate <- full_rate %>% fill(daily_rate, .direction="up")
  
  
  
  # FOR NSD/MSD  
  
  sub <- xx[2:15, ]
  centroid <- c(mean(coordinates(sub)[,1]), mean(coordinates(sub)[,2]))
  
  pj <- proj4::project(centroid, polar, inverse=T) #This transforms back the projected coordinates in geographic coordinates
  pj <- as.data.frame(pj)
  
  centroid <- c(pj$V1, pj$V2) #centroid in geographic coordinates (WGS84)
  
  xx$Loc_date <- format(as.Date(xx$Loc_date), "%Y-%m-%d") #Adding new time column w/ only the date no time 
  
  #Fonction suivante ne fonctionne pas si en SpatialPointsDataFrame...
  xx <- as.data.frame(xx)
  
  
  MeanInd <- as.data.frame(matrix(nrow = nrow(xx), ncol = 8))
  
  MeanInd <- xx  %>% #data frame used
    group_by(Loc_date) %>% #variable used to group the date for the mean
    summarise_at(vars(Longitude, Latitude), funs(mean(., na.rm=TRUE))) #specifiying for what we want to calculate a mean
  
  #Transforming projected to geographic since fct used to calculate distance use geographic
  
  matrix <- dplyr::select(MeanInd, Longitude, Latitude) #selecting projected coordinates and putting them in a matrix
  matrix <- as.matrix(matrix)
  
  matrixGeo <- proj4::project(matrix, polar, inverse = T) #Transforming projected to geographic
  matrixGeo <- as.data.frame(matrixGeo)
  
  #Adding the mean geographic coordinates to the MeanInd data frame
  MeanInd$long <- matrixGeo$V1
  MeanInd$lat <- matrixGeo$V2
  
  SP <- dplyr::select(MeanInd, long, lat) #Keeping only geographic coordinates
  SP <- SpatialPoints(SP, proj4string = wgs1984.proj)
  
  dist <- spDistsN1(SP, centroid, longlat = T) #Calculating distance from centroid (in km)
  
  dist <- as.data.frame(dist)
  
  MeanInd$dist <- dist$dist #Adding the distances to MeanInd 
  
  MeanInd$NSD <- ((MeanInd$dist)^2) #Calculating NSD from distance
  
  #Calculating MSD
  
  MeanInd$MSD <- rollmean(MeanInd$NSD, k = 5, fill = NA) #Here k is the size of the windows for which data are averaged 
  
  #Adding an ID and status collumn for the merging
  MeanInd$Fox_ID <- c(fox_ID)
  MeanInd$Loc_date <- as.POSIXct(strptime(xx$Loc_date,format="%Y-%m-%d", tz = "GMT"))
  
  
  #Create dataframe with all date
  xx$Loc_date <- as.POSIXct(strptime(xx$Loc_date,"%Y-%m-%d",tz = "GMT"))
  ts <- seq.POSIXt(xx$Loc_date[1], xx$Loc_date[nrow(xx)], by="day")
  df <- data.frame(Loc_date=ts)
  
  full_msd <- full_join(MeanInd, df)
  full_msd <- full_msd[order(full_msd[,"Loc_date"]),]
  full_msd$tracking_day <- seq.int(from=0,(nrow(full_msd)-1))
  
  #Fill missing rate value (giving the rate value that succeded the missing loc)
  full_msd <- full_msd %>% fill(MSD, .direction="up")
  
  
  table_id_disp <- as.data.frame(matrix(ncol = 15, nrow = nrow(full_rate)))
  
  colnames(table_id_disp) <- c("Fox_ID", "Platform_ID", "life_stage", "sex", "year","month", "loc_date", "tracking_day", "distance_Km", "daily_rate", "In0_Out1_Land", "longitude", "latitude", "NSD", "MSD")
  
  table_id_disp[,1] <- c(fox_ID)
  table_id_disp[,2] <- c(id_disp[[i]])
  table_id_disp[,3] <- c(life_stage)
  table_id_disp[,4] <- c(sex)
  table_id_disp[,5] <- as.factor(substr(full_rate$Loc_date,1,4))
  table_id_disp[,6] <- as.factor(substr(full_rate$Loc_date,6,7))
  
  
  letter.month <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  
  table_id_disp$month <- letter.month[table_id_disp$month]

  table_id_disp$month <- as.factor(table_id_disp$month)
  table_id_disp$month <- factor(table_id_disp$month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  
  
  
  table_id_disp[,7] <- full_rate$Loc_date
  table_id_disp[,8] <- full_rate$tracking_day
  table_id_disp[,9] <- full_rate$Distance_Km
  table_id_disp[,10] <-full_rate$daily_rate
  table_id_disp[,11] <-full_rate$In0_Out1_Land
  table_id_disp[,12] <-full_msd$long
  table_id_disp[,13] <-full_msd$lat
  table_id_disp[,14] <-full_msd$NSD
  table_id_disp[,15] <-full_msd$MSD
  
  
  
  table_daily_disp <- rbind(table_daily_disp, table_id_disp)
  
  #Start and end in tracking day of dispersal event
  start.disp <- table_id_disp$tracking_day[table_id_disp$loc_date==as.POSIXct(yy$start_disp,"%Y-%m-%d",tz = "GMT")]
  
  #if(is.na(yy$end_disp)){ 
    #end.disp <- table_id_disp$tracking_day[table_id_disp$Loc_date==as.POSIXct(ptt_id$end_tracking,"%Y-%m-%d",tz = "GMT")]
 # } else { 
 #   if(ptt_id$end_tracking=="") { 
      #end.disp <- table_id_disp$tracking_day[table_id_disp$Loc_date==as.POSIXct((xx$Loc_date[nrow(xx)]),"%Y-%m-%d",tz = "GMT")]
  #  } else { 
  #    end.disp <- table_id_disp$tracking_day[table_id_disp$Loc_date==as.POSIXct(yy$end_disp,"%Y-%m-%d",tz = "GMT")]
 #   }
  #} 
  
      end.disp <- table_id_disp$tracking_day[table_id_disp$loc_date==as.POSIXct(yy$end_disp,"%Y-%m-%d",tz = "GMT")]

  
  
  MSD.sansNA <- table_id_disp$MSD[!is.na(table_id_disp$MSD)]
  yaxis <- 150/max(MSD.sansNA)
  dist.sansNA <- table_id_disp$distance_Km[!is.na(table_id_disp$distance_Km)]
  dist.max <- round(max(dist.sansNA),0)
  pos.dist.max <- length(table_id_disp$tracking_day)*(1/8)

  
  gMerge <- ggplot(table_id_disp, aes(x=tracking_day))+
    geom_line(aes(y=rollmean(MSD, 5, na.pad=T)*yaxis), color="grey", size=2)+
    geom_line(aes(y=rollmean(daily_rate, 5, na.pad=T)), color="black", size=1)+
    scale_y_continuous(name="Daily rate (km/day)", sec.axis = sec_axis( ~./yaxis, name="MSD"), limits = c(0,150))+
    geom_vline(xintercept=c(start.disp), linetype="dotted", color="darkgreen", size=1)+
    geom_vline(xintercept=c(end.disp), linetype="dotted", color="darkred", size=1)+
    #annotate("text", x = pos.dist.max, y = 140, label = paste("Maximal distance =", dist.max, "km"),size=5)+
    xlab("Tracking days")+
    ggtitle(paste(fox_ID, id_disp[[i]], sep="_"))+
    guides(shape=TRUE)+
    theme_bw()+
    theme(axis.title.y.right=element_text(size=18,color="grey"),
          axis.title.y.left=element_text(size=18,color="black"),
          axis.title.x = element_text(size=18),
          legend.justification=c(0,0), legend.position=c(0,0))

  
  plot(gMerge)

}

#Export table
write.table(table_daily_disp,paste("/Table_daily_disp.txt"),col.names = T,row.names = F,sep = "\t")



#project database 

#Delete NAs (date added on missing days)
xxdisp <- subset(table_daily_disp, !is.na(longitude))

fox_disp_sp <- SpatialPointsDataFrame(xxdisp[, 12:13], data = xxdisp[, -c(12:13)], proj4string = wgs1984.proj)
fox_disp_sp <- spTransform(fox_disp_sp, polar)
str(fox_disp_sp)



#TABLE: parameters of dispersal events
#Create results matrix
disp <- as.data.frame(matrix(ncol = 27, nrow = length(id_disp)))
colnames(disp) <- c("Fox_ID", "Platform_ID", "life_stage", "sex", "start_tracking", "start_disp", "end_disp", "end_tracking", 
                    "disp_status", "end_status", "total_transfer_duration", "check_duration", "dist_tot", "dist_net", "dist_max","dist_final",
                    "linearity", "bearing", "bearing50", "daily_rate_emigration", "daily_rate_transfer", "daily_rate_immigration",
                    "daily_rate_transfer_land", "daily_rate_transfer_ice", "prop_land", "prop_ice", "complete")





for (i in 1 : length(id_disp)){
  
  xx <- fox_disp_sp[fox_disp_sp$Platform_ID == id_disp[[i]] ,]
  yy <- cols[cols$Platform_ID== id_disp[[i]] ,]
  #bear50 <- tab_bear50_sp[tab_bear50_sp$Platform_ID== id_disp[[i]] ,]

  
  #Duration of the transfer phase (complete or incomplete phase depending if establishment occurred)
  if(is.na(yy$end_disp)) { 
    da.end <- as.POSIXct(yy$end_tracking,"%Y-%m-%d",tz = "GMT")
  }else {
    da.end <- as.POSIXct(yy$end_disp,"%Y-%m-%d",tz = "GMT")
  }
  
  da.start <- as.POSIXct(yy$start_disp,"%Y-%m-%d", tz = "GMT")
  time_lag <- abs(difftime(da.start, da.end, units='days')) 
  disp[i,11] <- as.numeric(time_lag)
  
  

  #subset locs during transfer
  xxtrans <-xx[((substr(as.character(xx[,]$loc_date),1,10)) >= da.start & (substr(as.character(xx[,]$loc_date),1,10)) <= da.end),]
  #Check (number of locs during transfer)
  disp[i,12] <- nrow(xxtrans)
  
  #Fox, Platform, life stage and sex
  disp[i,1] <- levels(as.factor(as.character(xx$Fox_ID)))
  disp[i,2] <- levels(as.factor(as.character(xx[1,]$Platform_ID)))
  
  
  disp[i,3] <- levels(as.factor(as.character(yy$life_stage)))
  disp[i,4] <- levels(as.factor(as.character(yy$sex)))
  
  
  
  
  #tracking starting date
  disp[i,5] <- as.character(yy$start_tracking)
  
  #dispersal start/end date; establishment date for complete dispersal with an immigration phase
  disp[i,6] <- as.character(yy$start_disp)
  disp[i,7] <- as.character(yy$end_disp)
  
  #tracking ending date
  disp[i,8] <- as.character(yy$end_tracking)
  
  
  #long or short dispersal (short: did not >=80 km or >=30 days)
  disp[i,9] <- as.character(yy$disp_status)
  
  
  #ending of tracking status (A: still alive; B: collar removed; C: collar stopped; D: death of the animal)
  disp[i,10] <- as.character(yy$end_status)
  

  #Net distance (km) (eucledien distance between first and last locs of dispersal)
  #dist.ref 10, so dist.net is from end of buffer till last loc
  #dist.ref <- 5
#  dist.ref <- xxtrans[1,]$distance_Km
#  dist.final <- xxtrans[nrow(xxtrans),]$distance_Km
#  dist.net <- dist.final - dist.ref
#  disp[i,14] <- dist.net
  
#ERROR IN DIST MEASURE
  
  dist.ref <- 5
  dist.net <- (xxtrans[nrow(xxtrans),]$distance_Km) - dist.ref
  disp[i,14] <- dist.net
  
  
  
  
  #Max distance (km)
  disp[i,15] <- max(xx$distance_Km)
  
  #Final distance (km)
  disp[i, 16] <- (xx[nrow(xx),]$distance_Km) - dist.ref
  
  
 
  
  #Total distance (km)
  distance <- rep(NA,nrow(xxtrans@coords))
  ppp_points <- as(as(xxtrans, "SpatialPoints"), "ppp")  # convert point format
  
  for (j in 1 : length(distance)){
    Dist <- nncross(ppp_points[j,], ppp_points[j-1,])$dist    
    distance[j] <- Dist/1000   
  }
  
  distance[1] <-  0
  dist.tot <- sum(distance)
  disp[i,13] <- dist.tot
  
  #linearity
  disp[i,17] <- round((dist.net/dist.tot), digits = 4)
  
  
  #bearing 
  XXtrans<- proj4::project(xxtrans@coords, polar, inverse = T)
  #CHECK: a and f argument in bearing function
  bear <- bearing(XXtrans[1,], XXtrans[nrow(XXtrans),])
  disp[i,18] <- round((bear), digits = 3)
  
  
  #bearing 50km
  if(yy$Fox_ID == "RROB") { 
    disp[i,19] <- NA
  }else {
    row50km <- xx[xx$distance_Km >=50,]
    sub50<- row50km[1,]
    sub50 <- proj4::project(sub50@coords, polar, inverse = T)
    #CHECK: a and f argument in bearing function
    bear <- bearing(XXtrans[1,], sub50[1,])
    disp[i,19] <- round((bear), digits = 3)
  }
  
    
  
  

  
  #subset locs before transfer
  xxemi <-xx[((substr(as.character(xx[,]$loc_date),1,10)) >= xx$loc_date[1] & (substr(as.character(xx[,]$loc_date),1,10)) <= da.start),]
  
  #Mean daily rate before transfer
  disp[i,20] <- round((mean(xxemi$daily_rate)), digits = 4)
  
  
  #Mean daily rate during transfer
  disp[i,21] <- round((mean(xxtrans$daily_rate)), digits = 4)
  
  
  #subset locs after transfer
  
  if(is.na(yy$end_disp)){ 
    disp[i,22] <- NA
    
  } else { 
    if(is.na(yy$end_tracking)){ 
      disp[i,22] <- NA
      
    } else { 
    end.track <- as.POSIXct(yy$end_tracking,"%Y-%m-%d",tz = "GMT")
    xximmi <-xx[((substr(as.character(xx[,]$loc_date),1,10)) >= da.end & (substr(as.character(xx[,]$loc_date),1,10)) <= end.track),]

  #Mean daily rate after transfer
  disp[i,22] <- round((mean(xximmi$daily_rate)), digits = 4)
    }
  }
  
  #Proportion of locs on land and on sea ice during transfer
  xx_in <- subset(xxtrans, In0_Out1_Land=="0")
  prop_in <- round((nrow(xx_in)/nrow(xxtrans)), digits = 4)
  disp[i,25] <- prop_in
  
  xx_out <-  subset(xxtrans, In0_Out1_Land=="1")
  prop_out <-  round((nrow(xx_out)/nrow(xxtrans)), digits = 4)
  disp[i,26] <- prop_out
  
  #Mean daily rate on sea ice and on land during transfer
  
  if(prop_in == "0"){ 
    disp[i,23] <- NA
  } else { 
    daily_mean_in <- mean(xx_in$daily_rate)
    disp[i,23] <- round((daily_mean_in), digits = 4)
  }
 
  if(prop_out == "0"){ 
    disp[i,24] <- NA
  } else { 
    daily_mean_out <- mean(xx_out$daily_rate)
    disp[i,24] <- round((daily_mean_out), digits = 4)
  }

  if(yy$end_status == "immigrant"){ 
    disp[i,27] <- "yes"
  } else { 
    disp[i,27] <- "no"
  }
  
} 



#Export table
write.table(disp,paste("/Table_review_disp.txt"),col.names = T,row.names = F,sep = "\t")




              