#######################################
#### Classification of adult foxes ####
#######################################

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
setwd("")
fox <- read.delim("argos_filtrees_byday_2007_2021.txt", header = TRUE)

fox$Loc_date <- as.POSIXct(strptime(fox$Loc_date, format = "%Y-%m-%d %H:%M:%S", tz = "GMT"))

str(fox)

#Create objects for projections
wgs1984.proj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
UTMZ17 <- CRS("+proj=utm +zone=17 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
polar <- CRS("+proj=stere +lat_0=90 +lat_ts=71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")


#Remove the juveniles and the individuals with less than 14 locs of the database

fox <- fox[!fox[,"Fox_ID"] == "BBRM",] #Juvenile
fox <- fox[!fox[,"Fox_ID"] == "BBVR",] #Juvenile
fox <- fox[!fox[,"Fox_ID"] == "BJOJ",] #Juvenile
fox <- fox[!fox[,"Fox_ID"] == "BOBJ",] #Juvenile
fox <- fox[!fox[,"Fox_ID"] == "BORJ",] #Juvenile
fox <- fox[!fox[,"Fox_ID"] == "JOOB",] #Juvenile
fox <- fox[!fox[,"Fox_ID"] == "JOOJ",] #Juvenile
fox <- fox[!fox[,"Fox_ID"] == "JROJ",] #Juvenile
fox <- fox[!fox[,"Fox_ID"] == "MBBM",] #Juvenile
fox <- fox[!fox[,"Fox_ID"] == "MJVM",] #Juvenile
fox <- fox[!fox[,"Fox_ID"] == "OJBR",] #Juvenile
fox <- fox[!fox[,"Fox_ID"] == "OOBV",] #Juvenile
fox <- fox[!fox[,"Fox_ID"] == "OJOR",] #Juvenile
fox <- fox[!fox[,"Fox_ID"] == "OORO",] #Juvenile
fox <- fox[!fox[,"Fox_ID"] == "ORRB",] #Juvenile
fox <- fox[!fox[,"Fox_ID"] == "RBMV",] #Insuffisant number of locs
fox <- fox[!fox[,"Fox_ID"] == "RROB",] #Juvenile
fox <- fox[!fox[,"Fox_ID"] == "RVMR",] #Juvenile
fox <- fox[!fox[,"Fox_ID"] == "RVMM",] #Insuffisant number of locs
fox <- fox[!fox[,"Fox_ID"] == "VJRB",] #Juvenile
fox <- fox[!fox[,"Fox_ID"] == "VOBJ",] #Juvenile
fox <- fox[!fox[,"Fox_ID"] == "VOOB",] #Juvenile
fox <- fox[!fox[,"Fox_ID"] == "VVBO",] #Juvenile
fox <- fox[!fox[,"Fox_ID"] == "VVOB",] #Juvenile


fox <- mutate(fox, Fox_ID = factor(Fox_ID))
fox <- mutate(fox, Platform_ID = factor(Platform_ID))

fox_adu <- fox
str(fox_adu)

table(fox_adu$Platform_ID)
table(droplevels(fox_adu)$Platform_ID)


#Project the data
fox_adu_sp <- SpatialPointsDataFrame(fox_adu[, 3:4], data = fox_adu[, -c(3:4)], proj4string = wgs1984.proj)
fox_adu_sp <- spTransform(fox_adu_sp, polar)
#fox_adu_sp <- spTransform(fox_adu_sp, UTMZ17)


#Create a list with the Platform ID of the foxes
list_Platform_adu= sort(unique(as.character(fox_adu_sp$Platform_ID)))
list_Platform_adu



############################ Step 2 #############################
#### Determinate the presence or not of a stable home range #####
#################################################################


# Determinate if the individual has a stable home range based on the mean of the 14 firsts locs after capture

#Create a dataframe to stock data
table_HR <- data.frame(ID = list_Platform_adu, over = 0, prop = 0)


#Loop by Platform_ID
for (i in 1:length(list_Platform_adu)) {
    
  #Find the mean location of the 14 firsts locs after capture
  sub <- (fox_adu_sp[fox_adu_sp$Platform_ID == list_Platform_adu[[i]], ][2:15, ])
  cen <- c(mean(coordinates(sub)[,1]), mean(coordinates(sub)[,2]))
  cen_sp <- SpatialPoints(data.frame(cen[1], cen[2]), proj4string =polar)
  
  
  #Select the 14 firsts locs after capture
  sub_sp <- fox_adu_sp[fox_adu_sp$Platform_ID == list_Platform_adu[[i]], ][2:15, ]

  fox_ID <- droplevels(sub$Fox_ID[1])
  
  #shapefile(sub_sp, filename = paste(fox_ID, list_Platform_adu[[i]], "14locs", sep="_"), overwrite= TRUE)
  
  
  #Create a buffer of 5km ray around mean loc and see how many of the firsts 14 locs are in it
  tampon <- raster::buffer(cen_sp, width = 5000)
  table_HR[table_HR$ID == list_Platform_adu[[i]], "over"] <- sum(sp::over(sub_sp, tampon), na.rm = TRUE)
  
  #shapefile(tampon, filename = paste(fox_ID, list_Platform_adu[[i]], "FIRSTbuffer", sep="_"), overwrite= TRUE)

}


#Table with the proportion of the 14 firsts locs into the buffer
table_HR <- mutate(table_HR, prop = over/14)

#Graph
qplot(table_HR$prop,
      geom="histogram",
      main = "", 
      xlab = "Proportion of locs inside initial buffer",  
      ylab = "Number of foxes",
      fill=I("blue"), 
      col=I("darkblue"), 
      alpha=I(.4),
      bins = 10
)

#Separate individuals that have a stable home range (HR_stable) of those who have not (HR_instable)
HR_stable_adu <-  table_HR %>% dplyr::filter(prop>=0.5)
HR_instable_adu <-  table_HR %>% dplyr::filter(prop< 0.5)

list_ID_stable_adu = sort(unique(as.character(HR_stable_adu$ID)))
list_ID_stable_adu

list_ID_instable_adu = sort(unique(as.character(HR_instable_adu$ID)))
list_ID_instable_adu


###################### Step 3 #######################
####### Individuals with unstable home range ########
#####################################################

#Visual verification

list_ID_instable_adu
#BMBB (92256): Establish a territory
#BMRV (104309): Likely territory
#JRMR (113041): Territory or dead? 
#JVBJ (92259): Territory in C1 valley, very large
#OBOV (113051): Less likely territory
#OVOV (113058): Less likely territory (buffer in water, still floatter)
#VBJJ (104331): Likely territory?! Multiples excursions


#################### Step 4 ######################
######### Distances of locs from center ##########
##################################################


#Foxes with stable HR only are used here
list_ID_stable_adu

#Create a dataframe to stock data
tab <- data.frame()

#Loop by ID

for (i in 1 : length(list_ID_stable_adu)){
  
  #Subset locations 
  XY <- fox_adu_sp[fox_adu_sp@data$Platform_ID == list_ID_stable_adu[[i]],]
  
  #Mean of the 14 firsts locs
  sub <- (fox_adu_sp[fox_adu_sp$Platform_ID%in%list_ID_stable_adu[[i]], ][2:15, ])
  cen <- c(mean(coordinates(sub)[,1]), mean(coordinates(sub)[,2]))
  #cen_sp <- SpatialPoints(data.frame(cen[1], cen[2]), proj4string =polar)
  sub_first <- sub[1, ]
  
  center <- as.data.frame(matrix(nrow = 1, ncol = 6))
  colnames(center) <- c("Fox_ID","Platform_ID","Longitude","Latitude","Loc_date","Loc_quality")
  center[,1:2] <- sub_first@data[1:2]
  center[,5:6] <- NA
  center[,3:4] <- cen
  
  center <- SpatialPointsDataFrame(coords = center[,3:4], data = center[,c(1:2,5:6)], proj4string = UTMZ17)
  
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

  res <- mat

  
  #Order matrix by date
  ord <- order(res[,"Loc_date"])
  res <- res[ord,]

  #Stock into dataframe
  tab <- rbind(tab, res)
}

tab

################# Step 5 ####################
######### Extraction of excursions ##########
#############################################

#Buffer selection
# d = distance (km) from center
#Select distance
d <- 5

#Write 0 if the point is at less than distance d from center 
pos <- which(tab[,"Distance_Km"] < d)
if (length(pos) > 0) tab[pos,"Distance_Km"] <- 0

tab <- tab[order(tab[,"Fox_ID"],tab[,"Loc_date"]),]

setwd("/Bylot_Outline")

##In_Out of Bylot
#Create object for database
xy <- fox_adu_sp

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


#Exctract excursions: At least 1 points >5km of the center + point before and point after (complete path)
#All points of an excursion will have the same number on the last column

#Work with "tab" used before (ordered with date)
mat <- data.frame()
mat2 <- data.frame()

for (i in 1 : length(list_ID_stable_adu)){
  
  dat <- tab[tab[,"Platform_ID"] == list_ID_stable_adu[[i]],]
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

tad <- as.data.frame(matrix(ncol = 2,nrow = length(list_ID_stable_adu)))
colnames(tad) <- c("Fox_ID","Nb_Excursion")

for (i in 1 : length(list_ID_stable_adu)){
  
  x <- tot[tot[,"Platform_ID"] == list_ID_stable_adu[[i]],]
  tad[i,"Platform_ID"] <- list_ID_stable_adu[[i]]
  tad[i,"Nb_Excursion"] <- max(x[,"Excursion"])
  tad[i,"Fox_ID"] <- as.character(droplevels(x$Fox_ID[1]))
} 


#Number excursion's ID
mat[,"Excursion_fox"] = paste(mat[,"Fox_ID"], mat[,"Platform_ID"], mat[,"Excursion"],sep = "_")


tot <- tot[order(tot[,"Fox_ID"],tot[,"Loc_date"]),]
mat2 <- mat2[order(mat2[,"Fox_ID"],mat2[,"Loc_date"]),]
mat <- mat[order(mat[,"Fox_ID"],mat[,"Loc_date"]),]

#Export table
#setwd()
#write.table(tad, paste("/Summary_Excursions_adulte.txt",sep=""),col.names = T,row.names = F,sep = "\t")


#Final excursions table
mat_2 <- full_join(tab2, mat)
mat_2 <- drop_na(mat_2)
mat_2 <- mat_2[order(mat_2[,"Fox_ID"],mat_2[,"Loc_date"]),]

#mat_2 <- as.data.frame(cbind(mat_2[,1:2],mat_2[,6:7],mat_2[,3:4],mat_2[,9],mat_2[,11],mat_2[,8]))
mat_2 <- as.data.frame(cbind(mat_2[,1:2],mat_2[,5:6],mat_2[,3:4],mat_2[,8],mat_2[,10],mat_2[,7]))
colnames(mat_2) <- c("Fox_ID", "Platform_ID", "Longitude","Latitude","Loc_date","Loc_quality","Distance_Km","Excursion_fox","In0_Out1_Bylot")       
head(mat_2)
str(mat_2)


#Export table
write.table(mat_2,paste("/Table_Excursions_Adulte.txt"),col.names = T,row.names = F,sep = "\t")



#Extract vector of all excursions (of all foxes)
id_adu <- levels(as.factor(as.character(mat_2[,"Excursion_fox"])))

#Create results matrix
long_adu <- as.data.frame(matrix(ncol = 13,nrow = length(id_adu)))
colnames(long_adu) <- c("Fox_ID","Platform_ID","Excursion_ID","Date_Start","Date_End","Duration","Check", "Distance_Max","Distance_Final","Dispersion","In0_Out1_Bylot","Longitude_last", "Latitude_last")

#Explanation of columns:
#Date_start = Date of the first loc of the excursion/migration/dispersal event. (the loc before leaving HR)
#Date_End = Date of the last loc of the excursion/migration/dispersal event. 
#Check = number of locs of the path (exclude locs into the HR), so we can see if there are locs missing in the path
#Duration = duration in days of the path (should normally match the "Check" column)
#Distance_Max = maximal distance (km) of the path from center
#Distace_Final = distance (km) of the last loc of the path
#Dispersion = If there is a return to the HR after the excursion(0), else classed into dispersal event(1)
#In0_Out1_Bylot = If the excursion/dispersion finish on Bylot(0) or outside Bylot(1)
#Longitude_last & Latitude_last = coordinates of the last loc of the excursion


for (i in 1 : length(id_adu)){
  
  xx <- mat_2[mat_2[,"Excursion_fox"] == id_adu[i],]
  
  #Fox, Platform and Excursion ID
  long_adu[i,1] <- levels(as.factor(as.character(xx$Fox_ID)))
  long_adu[i,2] <- levels(as.factor(as.character(xx[1,]$Platform_ID)))
  long_adu[i,3] <- levels(as.factor(as.character(xx$Excursion_fox)))
  
  #Begining and ending of the excursion/dispersion
    #Time lag between first and last loc if excursion
    long_adu[i,4] <- as.character(xx[1,]$Loc_date)
    long_adu[i,5] <- as.character(xx[nrow(xx),]$Loc_date)


  #Duration of the excursion
  if (xx[nrow(xx),]$Distance_Km == "0") { 
    #Time lag between first and last loc if excursion
    da.start <- as.POSIXct(strptime(as.character(xx[1,]$Loc_date),"%Y-%m-%d", tz = "GMT"))
    da.end <- as.POSIXct(strptime(as.character(xx[nrow(xx),]$Loc_date),"%Y-%m-%d",tz = "GMT"))
    time_lag <- (abs(difftime(da.start, da.end, units='days'))-1) #Taking of the day before departure into the HR
    long_adu[i,6] <- as.numeric(time_lag)
  } else {
    #Time lag between first and last loc if dispersion
    da.start <- as.POSIXct(strptime(as.character(xx[1,]$Loc_date),"%Y-%m-%d", tz = "GMT"))
    da.end <- as.POSIXct(strptime(as.character(xx[nrow(xx),]$Loc_date),"%Y-%m-%d",tz = "GMT"))
    time_lag <- abs(difftime(da.start, da.end, units='days')) #Taking of the day before departure into the HR
    long_adu[i,6] <- as.numeric(time_lag)
  }
  
  #Number of locs for each excursion
  if (xx[nrow(xx),]$Distance_Km == "0") { #excursion (last loc is in HR Distance = 0)
    if ((nrow(xx)-2) > 1) {
      long_adu[i,7] <- nrow(xx)-2 #nb of locs without locs of departure and arrival inside HR
    } else {
      long_adu[i,7] <- 1
    }
  } else { #dispersion (last loc isn't in HR Distance different of 0)
    long_adu[i,7] <- nrow(xx)-1 #nb of locs without loc of departure inside HR
  }
  
  #Max distance of the excursion
  long_adu[i,8] <- max(xx$Distance_Km-5)
  
  #Final distance (to use for dispersion)
  if (xx[nrow(xx),]$Distance_Km == "0") { 
    long_adu[i,9] <- (xx[nrow(xx)-1,]$Distance_Km)
  } else {
    long_adu[i,9] <- (xx[nrow(xx),]$Distance_Km)
  }
  
  #Separate excursion/migration of dispersion (0 if excursion, return into the HR; 1 if dispersion, did'nt return into HR)
  if (xx[nrow(xx),]$Distance_Km == "0") { 
    long_adu[i,10] <- 0
  } else {
    long_adu[i,10] <- 1
  }
  
  #Determinate if the excursion or dispersion finish on Bylot
  if (xx[nrow(xx),]$In0_Out1_Bylot == "0") { 
    long_adu[i,11] <- 0
  } else {
    long_adu[i,11] <- 1
  }

  #Add columns for last longitude and latitude data
  if (xx[nrow(xx),]$Distance_Km == "0") { 
    long_adu[i,12] <- (xx[nrow(xx)-1,]$Longitude)
  } else {
    long_adu[i,12] <- (xx[nrow(xx),]$Longitude)
  } 
  
  if (xx[nrow(xx),]$Distance_Km == "0") { 
    long_adu[i,13] <- (xx[nrow(xx)-1,]$Latitude)
  } else {
    long_adu[i,13] <- (xx[nrow(xx),]$Latitude)
  }  

}

long_adu


########################## Step 6 #############################
#### Classification of Residency, Migration and Dispersion ####
###############################################################


dist_max_80 <- long_adu %>% dplyr::filter(Distance_Max>=80) %>% filter(Distance_Final <80) %>% filter(Dispersion==1)
dist_max_80$Excursion_ID


# LONG EXCURSION EVENTS
#Excursions of more than 35 days AND/OR more than 80 km from HR center

long_excursions_adu <- long_adu %>% dplyr::filter(Distance_Final>=80) 
long_excursions_adu$Excursion_ID


# DISPERSAL EVENTS
#Long excursion with no return in HR

long_disp_adu <- long_adu %>% dplyr::filter(Distance_Max>=80) %>% filter(Dispersion==1)
long_disp_adu$Excursion_ID


# MIGRATORY EVENTS
#Long excursion with return in HR

long_mig_adu <- long_adu %>% dplyr::filter(Distance_Max>=80) %>% filter(Dispersion==0)
long_mig_adu$Excursion_ID

