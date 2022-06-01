###############################################
##### Classification Juvenile Dispersal #####
###############################################


# Based on Sandra Lai's script
# Modified by Richard Gravel

# Download libraries
library(readxl)
library(rgeos)
library(sp)
library(proj4)
library(dplyr)
library(spdplyr)
library(raster)
library(rgdal)
library(ggplot2)
library(adehabitatLT)
library(tidyr)
library(sf)
library(geosphere)
library(maptools)
library(spatstat)

################## Step 1 ###################
############## Database check ###############
#############################################

#Download database (one loc by day)
setwd( "")
fox <- read.delim("argos_filtrees_byday_2007_2021.txt", header = TRUE)
str(fox)


#Create objects for projections
wgs1984.proj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
UTMZ17 <- CRS("+proj=utm +zone=17 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
polar <- CRS("+proj=stere +lat_0=90 +lat_ts=71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

# Select individuals
fox_juv <- subset(fox, Fox_ID=="BBRM" | Fox_ID=="BBVR" | Fox_ID=="BOBJ" | Fox_ID=="JROJ" | Fox_ID=="RVMR" | Fox_ID=="VJRB" | Fox_ID=="VOBJ" | Fox_ID=="VOOB" | Fox_ID=="VVOB" | Fox_ID=="BJOJ" | Fox_ID=="BORJ" | Fox_ID=="JOOB" | Fox_ID=="OOBV" | Fox_ID=="OORO" | Fox_ID=="JOOJ"| Fox_ID=="OJBR" | Fox_ID=="OJOR" | Fox_ID=="ORRB" | Fox_ID=="RROB" | Fox_ID=="VVBO")

#Juveniles not included:
#MBBM - 7 locs only (loc_quality too bad, all of the 7 locs were left out)
#MJVM - lost collar

fox_juv <- fox_juv[order(fox_juv[,"Fox_ID"],fox_juv[,"Loc_date"]),]

fox_juv <- mutate(fox_juv, Fox_ID = factor(Fox_ID))
table(droplevels(fox_juv)$Fox_ID)
str(fox_juv)

#Project the data
fox_juv_sp <- SpatialPointsDataFrame(fox_juv[, 3:4], data = fox_juv[, -c(3:4)], proj4string = wgs1984.proj)
fox_juv_sp <- spTransform(fox_juv_sp, polar)
fox_juv_sp

#Create a list with the ID of the foxes
list_ID_juv = sort(unique(as.character(fox_juv_sp$Fox_ID)))
list_ID_juv

list_Platform_juv= sort(unique(as.character(fox_juv_sp$Platform_ID)))
list_Platform_juv


#################### Step 2 ######################
######### Distances of locs from center ##########
##################################################

#Buffer around the first loc. Juveniles are captured on the natal den.

#Create a dataframe to stock data
tab <- data.frame()

#Loop by ID
for (i in 1 : length(list_ID_juv)){
  
  XY <- fox_juv_sp[fox_juv_sp@data$Fox_ID == list_ID_juv[[i]],]
  
  
  #Center
  center <- (fox_juv_sp[fox_juv_sp$Fox_ID%in%list_ID_juv[[i]], ][1, ])
  
  
  # Loc's distance form center 
  #===================================
  
  # Calculation of minimal Eucledean distances from center
  
  distance <- rep(NA,nrow(XY@coords))
  ppp_points <- as(as(XY, "SpatialPoints"), "ppp")  # convert point format
  psp_center <- as.ppp.SpatialPointsDataFrame(center) # convert center in ppp format
  
  for (j in 1 : length(distance)){
    
    # Calcul de toutes les distances du point j avec le centre
    Dist <- nncross(ppp_points[j,], psp_center)$dist    
    distance[j] <- Dist        
  }
  
  #Create a matrix of distances for one individual
  #Matrix with distances of points from center
  mat <- as.data.frame(matrix(nrow = length(XY), ncol = 7))
  colnames(mat) <- c("Fox_ID","Platform_ID","Longitude","Latitude","Loc_date","Loc_quality","Distance_Km")
  mat[,1:2] <- XY@data[1:2]
  mat[,5:6] <- XY@data[3:4]
  mat[,3:4] <- XY@coords
  mat[,7] <- distance/1000
  
  #Matrix with center indicated dist=0
  mat2 <- as.data.frame(matrix(nrow = length(center), ncol = 7))
  colnames(mat2) <- c("Fox_ID","Platform_ID", "Longitude","Latitude","Loc_date","Loc_quality", "Distance_Km")
  mat2[,1:2] <- center@data[1:2]
  mat2[,5:6] <- center@data[3:4]
  mat2[,3:4] <- center@coords
  mat2[,7] <- 0
  
  #Merge the 2 matrices
  #?!?!?!??!?!?!?!?!??!?!
  ##res <- rbind(mat,mat2)
  res <- mat
  
  
  #Order matrix by date
  
  ord <- order(res[,"Loc_date"])
  res <- res[ord,]
  
  #Stock into dataframe
  tab <- rbind(tab, res)
}


################# Step 3 ####################
######### Extraction of excursions ##########
#############################################

#Buffer selection
# d = distance (km) from center
#Select distance
###Dans code precedent, j'ai utilise 20 km. alors que j'aurais du selectionner 10?!
d <- 5

#Write 0 if the point is at less than distance d from center 
pos <- which(tab[,"Distance_Km"] < d)
if (length(pos) > 0) tab[pos,"Distance_Km"] <- 0

tab <- tab[order(tab[,"Fox_ID"],tab[,"Loc_date"]),]

# Set directiory
setwd("/Bylot_Outline")
##In_Out of Bylot

#Create object for database
xy <- fox_juv_sp

#Import Bylot outline
Bylot <- rgdal::readOGR("Bylot_outline_polar.shp")

#Polygon into SpatialPolygons format
bylot_pol <- Bylot[,1]@polygons
bylot <- SpatialPolygons(bylot_pol, proj4string = polar)

# Clip Bylot on points for all foxes
clipped <- over(xy, bylot)

# Points IN Bylot
xy_in <- xy[which(!is.na(clipped)),]

# Points OUT Bylot
xy_out <- xy[which(is.na(clipped)),]

dat <- as.data.frame(matrix(nrow = length(xy_out), ncol = 1))
In0_Out1Bylot <- rep("1", nrow(xy_out))
comp <- as.data.frame(xy_out)
yx2 <- cbind(comp,In0_Out1Bylot)

dat <- as.data.frame(matrix(nrow = length(xy_in), ncol = 1))
In0_Out1Bylot <- rep("0", nrow(xy_in))
comp <- as.data.frame(xy_in)
xy_in2 <- cbind(comp,In0_Out1Bylot)

tab2 <- rbind(yx2, xy_in2)
ord <- order(tab2[,"Fox_ID"])
tab2 <- tab2[ord,]


#Exctract excursions: At least 1 points >20km of the center + point before and point after (complete path)
#All points of an excursion will have the same number on the last column

#Work with "tab" used before (ordered with date)
mat <- data.frame()
mat2 <- data.frame()

for (i in 1 : length(list_ID_juv)){
  
  dat <- tab[tab[,"Fox_ID"] == list_ID_juv[[i]],]
  rownames(dat) <- NULL
  dat[,"Excursion"] <- 0
  k <- 1
  for (j in 1 : nrow(dat)){
    
    if (j < nrow(dat)){
      if (dat[j,"Distance_Km"] != 0){
        dat[j,"Excursion"] <- k
        if (dat[j+1,"Distance_Km"] == 0){
          k <- k + 1
        }
      }
    }else{
      if (dat[j,"Distance_Km"] != 0){
        dat[j,"Excursion"] <- k
      }
    }
  }
  
  k <- max(dat[,"Excursion"])
  if (k > 0) {
    for (w in 1 : k){
      
      x1 <- x2 <- NULL
      pos <- which(dat[,"Excursion"] == w)
      
      if (pos[1] != 1) x1 <- pos[1]-1
      if (pos[length(pos)] != nrow(dat)) x2 <- pos[length(pos)]+1
      
      if (!is.null(x1) && !is.null(x2)) pos <- c(x1,pos,x2)
      if (is.null(x1) && !is.null(x2)) pos <- c(pos,x2)
      if (!is.null(x1) && is.null(x2)) pos <- c(x1,pos)
      
      y <- dat[pos,]
      y[,"Excursion"] <- w
      mat <- rbind(mat,y)
      rownames(mat) <- NULL
    }  
  }else{
    mat2 <- rbind(mat2,dat)
  } 
}



# Number of excursion by fox
tot <- rbind(mat, mat2)

tad <- as.data.frame(matrix(ncol = 2,nrow = length(list_ID_juv)))
colnames(tad) <- c("Fox_ID","Nb_Excursion")

for (i in 1 : length(list_ID_juv)){
  
  x <- tot[tot[,"Fox_ID"] == list_ID_juv[[i]],]
  tad[i,"Fox_ID"] <- list_ID_juv[[i]]
  tad[i,"Nb_Excursion"] <- max(x[,"Excursion"])
}

tad_juv <- tad

#Number excursion's ID
mat[,"Excursion_fox"] = paste(mat[,"Fox_ID"], mat[,"Excursion"],sep = "_")

#Export tables
#setwd("/Tables_Excursions")
#write.table(tad_juv, paste("/Summary_Excursions_Juv.txt",sep=""),col.names = T,row.names = F,sep = "\t")


#Final excursions table
mat_2 <- full_join(tab2, mat)
mat_2 <- drop_na(mat_2)
mat_2 <- mat_2[order(mat_2[,"Fox_ID"],mat_2[,"Loc_date"]),]

mat_2 <- as.data.frame(cbind(mat_2[,1:2],mat_2[,5:6],mat_2[,3:4],mat_2[,8],mat_2[,10],mat_2[,7]))
colnames(mat_2) <- c("Fox_ID", "Platform_ID", "Longitude","Latitude","Loc_date","Loc_quality","Distance_Km","Excursion_fox","In0_Out1_Bylot")       
head(mat_2)
str(mat_2)

mat_2_juv <- mat_2

#Export table
write.table(mat_2_juv,paste("/Users/richardgravel/Documents/Manuscrit/data/Table_Excursions_Juv.txt"),col.names = T,row.names = F,sep = "\t")


#Excursion Table

#Extract vector of all excursions (of all foxes)
id_juv <- levels(as.factor(as.character(mat_2[,"Excursion_fox"])))

#Create results matrix
long_juv <- as.data.frame(matrix(ncol = 13,nrow = length(id_juv)))
colnames(long_juv) <- c("Fox_ID","Platform_ID","Excursion_ID","Date_Start","Date_End","Duration","Check", "Distance_Max","Distance_Final","Dispersion","In0_Out1_Bylot","Longitude_last", "Latitude_last")


for (i in 1 : length(id_juv)){
  
  xx <- mat_2[mat_2[,"Excursion_fox"] == id_juv[i],]
  
  #Fox, Platform and Excursion ID
  long_juv[i,1] <- levels(as.factor(as.character(xx$Fox_ID)))
  long_juv[i,2] <- levels(as.factor(as.character(xx[1,]$Platform_ID)))
  long_juv[i,3] <- levels(as.factor(as.character(xx$Excursion_fox)))
  
  
  #Begining and ending of the excursion/dispersion
  
  #Time lag between first and last loc if excursion
  long_juv[i,4] <- as.character(xx[1,]$Loc_date)
  long_juv[i,5] <- as.character(xx[nrow(xx),]$Loc_date)
  
  
  #Duration of the excursion
  if (xx[nrow(xx),]$Distance_Km == "0") { 
    #Time lag between first and last loc if excursion
    da.start <- as.POSIXct(strptime(as.character(xx[1,]$Loc_date),"%Y-%m-%d", tz = "GMT"))
    da.end <- as.POSIXct(strptime(as.character(xx[nrow(xx),]$Loc_date),"%Y-%m-%d",tz = "GMT"))
    time_lag <- (abs(difftime(da.start, da.end, units='days'))-1) #Taking of the day before departure into the HR
    long_juv[i,6] <- as.numeric(time_lag)
  } else {
    #Time lag between first and last loc if dispersion
    da.start <- as.POSIXct(strptime(as.character(xx[1,]$Loc_date),"%Y-%m-%d", tz = "GMT"))
    da.end <- as.POSIXct(strptime(as.character(xx[nrow(xx),]$Loc_date),"%Y-%m-%d",tz = "GMT"))
    time_lag <- abs(difftime(da.start, da.end, units='days')) #Taking of the day before departure into the HR
    long_juv[i,6] <- as.numeric(time_lag)
  }
  
  #Number of locs for each excursion
  if (xx[nrow(xx),]$Distance_Km == "0") { #excursion (last loc is in HR Distance = 0)
    if ((nrow(xx)-2) > 1) {
      long_juv[i,7] <- nrow(xx)-2 #nb of locs without locs of departure and arrival inside HR
    } else {
      long_juv[i,7] <- 1
    }
  } else { #dispersion (last loc isn't in HR Distance different of 0)
    long_juv[i,7] <- nrow(xx)-1 #nb of locs without loc of departure inside HR
  }
  
  #Max distance of the excursion
  long_juv[i,8] <- max(xx$Distance_Km)
  
  #Final distance of the excursion/dispersal
  #Pour une excursion... ca n'a pas vraiment de valeur. faudrait avoir cette colonne pour dispersion seulement
  
  if (xx[nrow(xx),]$Distance_Km == "0") { 
    long_juv[i,9] <- (xx[nrow(xx)-1,]$Distance_Km)
  } else {
    long_juv[i,9] <- (xx[nrow(xx),]$Distance_Km)
  }
  
  
  #Separate excursion/migration of dispersion (0 if excursion, return into the HR; 1 if dispersion, did'nt return into HR)
  if (xx[nrow(xx),]$Distance_Km == "0") { 
    long_juv[i,10] <- 0
  } else {
    long_juv[i,10] <- 1
  }
  
  #Determinate if the excursion or dispersion finish on Bylot
  if (xx[nrow(xx),]$In0_Out1_Bylot == "0") { 
    long_juv[i,11] <- 0
  } else {
    long_juv[i,11] <- 1
  }
  
  
  #Add columns for last longitude and latitude data
  if (xx[nrow(xx),]$Distance_Km == "0") { 
    long_juv[i,12] <- (xx[nrow(xx)-1,]$Longitude)
  } else {
    long_juv[i,12] <- (xx[nrow(xx),]$Longitude)
  } 
  
  if (xx[nrow(xx),]$Distance_Km == "0") { 
    long_juv[i,13] <- (xx[nrow(xx)-1,]$Latitude)
  } else {
    long_juv[i,13] <- (xx[nrow(xx),]$Latitude)
  }  
  
}

long_juv


#Join BD with tracking info
#tracking_info <- argos_periode_suivi
#colnames(tracking_info) <- c("noid", "Platform_ID","ptt_cls","date_debut_suivi","date_fin_suivi","Tracking_Status","notes")

#long_juv <- merge(x = long_juv, y = tracking_info, by = "Platform_ID")
#long_juv$noid <- long_juv$ptt_cls <- long_juv$date_debut_suivi <- long_juv$date_fin_suivi <-  NULL


#Order table by Excursion ID
long_juv <- long_juv[order(long_juv[,"Excursion_ID"], long_juv[,"Date_Start"]),]


########################## Step 4 #############################
#### Classification of Residency, Migration and Dispersion ####
###############################################################

# LONG EXCURSION EVENTS

#Excursions of more than 80 km from HR center

long_excursions_juv <- long_juv %>% dplyr::filter(Distance_Max>=80) 
long_excursions_juv$Excursion_ID
# [1] "BBVR_1" "BOBJ_2" "BORJ_2" "JOOB_9" "JOOJ_1" "OJBR_1" "OJOR_6" "OOBV_4" "OOBV_6" "OORO_1"
# [11] "VOBJ_5" "VOOB_4" "VVOB_3"



# DISPERSAL EVENTS
#Excursion with no return in HR
long_dispersal_juv <- long_juv %>% dplyr::filter(Distance_Final>=80) %>% filter(Dispersion==1)
long_dispersal_juv$Excursion_ID
#   "BBVR_1" "BOBJ_2" "BORJ_2" "JOOB_9" "JOOJ_1" "OJBR_1" "OJOR_6" "OORO_1" "VOOB_4" "VVOB_3"


# MIGRATORY EVENTS
#Long excursion with return in HR

long_mig_juv <- long_juv %>% dplyr::filter(Distance_Max>=80) %>% filter(Dispersion==0)
long_mig_juv$Excursion_ID
#[1] "OOBV_4" "OOBV_6" "VOBJ_5"

#notes:
#VOBJ: dispersal afterward?



