################################################################################
# swcFilter.v02.R (edited from July to August 2014 
################################################################################
#
# USAGE: This code removes the failed SWC measurements produced by CS610 soil
# 	 water content probes located at PJ sites from 2009 to 2013. 
#
#     These data have been characterized
#	 as having a consistent 'electrical noise' signal which results in a 
#	 seemingly random downward spike in the data at intermittent intervals
#	 but apparently consistent magnitudes. Given the desire to preserve the
# 	 perceived 'real' positive spikes in the signal due to precipitation,
#	 but not negatively affect the mean signal during dry, noisey periods,
#	 two stages filter is applied:
#
#  First Stage: a wavelet filter may produce better results than a median, variance,
#	 or mean moving average window (time domain) based approach.It generates OUTPUT 1.
#
#  Second stage:remaining errors in SWC measurements are filtered by using a quality (QC) 
#  indicator that presents 3 different forms from less to more restrictive. This second 
#  stage precises of manual decissions to be made in order to remove only problematic data.
#  It generates OUTPUT 2.

#                  QC_dif (orange):  SWC(i)-SWC(i-1) *100/ SWC(i-1)       % of variation     (+ when SWC is increasing  and - when SWC is decreasing )
#
#                  QC_2 (green): If it rained the last half a day and QC_diff is positive QC_2=0, otherwise QC_2= QC_diff

#                  QC_3 (purple):
#                      If  rain =0 (from last day) and QC_dif>0 ? QC_3= 1
#                      If  SWC is decreasing more rapidly than the maximum REALISTIC SWC decreased (QC_2 >DRY MAX) ?  QC_3= 1
#                      Otherwise ?QC_3= 0

#                  QC_4 (red): standard deviation of QC_3 during a time window of ¾ of a day

#                  QC_5(brown): if QC_4> 0.5 max(QC_4)?QC_5=1, Oherwise QC_5=0

#                  QC_6(pink): if QC_4> 0 ?QC_6=1, Oherwise QC_6=0

#   This script also includes one last section to estimate the daily means of SWC and also to gapfill the daily gaps when there
#   was not rain during the gap and neither for the first and last day when SWC is available


#    6 OUTPUT FILES:
#        > OUTPUT1_name :  "SWC_site_year_Median_filtered_30min"
#        > OUTPUT2_name: "SWC_site_year_REFINED_30min"
#        > OUTPUT3_name: "SWC_site_year_REFINED_Daily"
#        > OUTPUT4_name: "SWC_site_year_Refined&Gapfilled_Daily"
#        > OUTPUT5_name: "SWC_site_year_Refined&Gapfilled_Daily_with_means"
#         GAPS_SWC_PJC_20010"
# 
#    4 pdf DIAGNOSTIC PLOTS:
#        > PLOT1_name: "PlotSWC_site_year_Median_filtered_30min"
#        > PLOT2_name: "PlotSWC_site_year_REFINED_30min"
#        > PLOT3_name: "PlotSWC_site_year_REFINED_Daily"
#        > PLOT4_name: "PlotSWC_site_year_Refined&Gapfilled_Daily"

        
# For inquiries about this script's use or function, contact the authors at:
#	Dan Krofcheck      krofcheck@gmail.com
# Laura Morillas     laura.morillas.gonzalez@gmail.com
#
#################################################################################
#################################################################################

#For naming OUTPUT files:

SITES<-c("PJC","PJG")

SITE.CODE<-2              #WRITE HERE THE CODE (1 for PJC, 2 for PJG) OF THE SITE TO PROCESS
YEAR<-"2010"              #WRITE HERE THE YEAR OF YEARS OF DATA TO PROCESS

SITE<-SITES[SITE.CODE]
LEVEL1<-"Median_filtered"
LEVEL2<-"REFINED"
LEVEL3<-"Refined&Gapfilled"
TSCALE1<-"30min"
TSCALE2<-"Daily"

OUTPUT1_name<-paste("SWC",SITE,YEAR,LEVEL1,TSCALE1, sep="_")
OUTPUT2_name<-paste("SWC",SITE,YEAR,LEVEL2,TSCALE1, sep="_")
OUTPUT3_name<-paste("SWC",SITE,YEAR,LEVEL2,TSCALE2, sep="_")
OUTPUT4_name<-paste("SWC",SITE,YEAR,LEVEL3,TSCALE2, sep="_")
OUTPUT5_name<-paste("SWC",SITE,YEAR,LEVEL3,TSCALE2,"with_means", sep="_")

PLOT1_name<-paste("PlotSWC",SITE,YEAR,LEVEL1,TSCALE1, sep="_")
PLOT2_name<-paste("PlotSWC",SITE,YEAR,LEVEL2,TSCALE1, sep="_")
PLOT3_name<-paste("PlotSWC",SITE,YEAR,LEVEL2,TSCALE2, sep="_")
PLOT4_name<-paste("PlotSWC",SITE,YEAR,LEVEL3,TSCALE2, sep="_")

GAP_name<-paste("GAPS_SWC",SITE,YEAR,sep="_")

#################################
##
##Location, functions and Packages
##
#################################

setwd("C:/Users/Laura Morillas/Documents/R/SWC filtering")

### Depends package Reshape

#install.packages("ggplot2")
#install.packages("reshape")
#install.packages("robfilter")

library(ggplot2)
library(reshape)
library(zoo)
library(robfilter)

myfx<- function (x){ sum(!is.na(x))} #function to account the number of values that are not NA

################################################
#IMPORTING DATA
################################################

### input file will be a CSV file containing raw SWC data 

#PJC 2009:
rawData <- read.table("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/PJC_soil_complete/PJ_2009_soil_complete.dat", sep = "\t", header = TRUE)
fluxAll <- read.table("C:/Users/Laura Morillas/Documents/data/PJ_control/Flux_all files/PJ_FLUX_all_2009.txt",sep = "\t", header = TRUE)
fluxAll<-fluxAll[,c(1,3:165)]

#PJC 2010:
rawData <- read.table("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/PJC_soil_complete/PJ_2010_soil_complete.dat", sep = "\t", header = TRUE)
fluxAll <- read.table("C:/Users/Laura Morillas/Documents/data/PJ_control/Flux_all files/PJ_FLUX_all_2010.txt",sep = "\t", header = TRUE)
fluxAll<-fluxAll[1:nrow(fluxAll)-1,2:ncol(fluxAll)]


#PJC 2011:
rawData <- read.table("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/PJC_soil_complete/PJ_2011_soil_complete_07242014.txt", sep = ",", header = TRUE)
rawData<-PJC_rawData[2:17520,]            #only for 2011 data because:all soil_Year_complete datasets begin doy 1 00:30 and finish doy 365 or 366 at 23:30
fluxAll  <- read.csv("C:/Users/Laura Morillas/Documents/data/PJ_control/Flux_all files/PJ_FLUX_all_2011_Dan.csv",header=T,sep = ",",dec=".")
fluxAll<-PJC_fluxAll[2:nrow(fluxAll),]


#PJC 2012:
rawData <- read.table("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/PJC_soil_complete/PJ_2012_soil_complete_07242014.txt", sep = ",", header = TRUE)
rawData<-rawData[2:17568,]            #only for 2011 data because:all soil_Year_complete datasets begin doy 1 00:30 and finish doy 365 or 366 at 23:30

fluxAll<- read.table("C:/Users/Laura Morillas/Documents/data/PJ_control/Flux_all files/PJ_FLUX_all_2012.txt",sep ="\t", header = TRUE)
fluxAll<- fluxAll[1:17567,1:152]


#PJC 2013:
rawData <- read.table("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/PJC_soil_complete/PJ_2013_soil_complete - Randy_Lauracompleted.txt", sep = "\t", header = TRUE)

fluxAll  <- read.csv("C:/Users/Laura Morillas/Documents/data/PJ_control/Flux_all files/PJ_FLUX_all_2013.txt",header=T,sep = "\t",dec=".")
fluxAll<-fluxAll[2:nrow(fluxAll),1:152]



#PJG 2011:
FLUX_SOIL<-read.table("C:/Users/Laura Morillas/Documents/data/PJ_guirdled/soil/From bandage script/PJ_girdle_2011_corrected_tstamps_23x_fluxall.dat", sep = "\t", header = TRUE) #sep = "\t" means separated by Tab
rawData<-FLUX_SOIL[1:(nrow(FLUX_SOIL)-1),162:ncol(FLUX_SOIL)]   #to have only data coming from CR23X (FROM COLUMN 162 OF THE ORIGINAL FILE) including data doy 1 00:30 and finish doy 365 or 366 at 23:30 
                                                             # (original soil complete time columns) 
            
                                                           
fluxAll <- read.table("C:/Users/Laura Morillas/Documents/data/PJ_guirdled/Flux_all files/PJ_girdle_FLUX_all_2011.txt",sep = "\t", header = TRUE)
fluxAll<-fluxAll[1:nrow(fluxAll)-1,]   # To have Fux all file from doy 1 00:30 and finish doy 365 or 366 at 23:30 
#We need some PJC 2011 DATA TO SUPPORT pjg 2011 DATA.....BECAUSE TIMESTAMP AT PJG IS NOT OK
PJC_rawData <- read.table("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/PJC_soil_complete/PJ_2011_soil_complete_07242014.txt", sep = ",", header = TRUE)
PJC_rawData<-PJC_rawData[2:17520,]            #only for 2011 data because:all soil_Year_complete datasets begin doy 1 00:30 and finish doy 365 or 366 at 23:30
PJC_fluxAll  <- read.csv("C:/Users/Laura Morillas/Documents/data/PJ_control/Flux_all files/PJ_FLUX_all_2011_Dan.csv",header=T,sep = ",",dec=".")
PJC_fluxAll<-PJC_fluxAll[2:nrow(PJC_fluxAll),]


#PJG 2009:
rawData <- read.table("C:/Users/Laura Morillas/Documents/data/PJ_guirdled/soil/PJG_soil_complete/PJ_girdle_2009_soil_complete.dat", sep = "\t", header = TRUE)
fluxAll <- read.table("C:/Users/Laura Morillas/Documents/data/PJ_guirdled/Flux_all files/PJ_girdle_FLUX_all_2009.txt",sep ="\t", header = TRUE)
fluxAll<-fluxAll[,2:159]
rawData[,1]<-PJC_rawData[,1]



#PJG 2010:
FLUX_SOIL<-read.table("C:/Users/Laura Morillas/Documents/data/PJ_guirdled/soil/From bandage script/PJ_girdle_2010_corrected_tstamps_23x_fluxall.dat", sep = "\t", header = TRUE) #sep = "\t" means separated by Tab
rawData<-FLUX_SOIL[1:(nrow(FLUX_SOIL)-1),162:ncol(FLUX_SOIL)]   #to have only data coming from CR23X (FROM COLUMN 162 OF THE ORIGINAL FILE) including data doy 1 00:30 and finish doy 365 or 366 at 23:30 
#timestampDF<- rawData[,2:4]

fluxAll <- read.table("C:/Users/Laura Morillas/Documents/data/PJ_guirdled/Flux_all files/PJ_girdle_FLUX_all_2010.txt",sep = "\t", header = TRUE)
fluxAll<-fluxAll[1:nrow(fluxAll)-1,]   # To have Fux all file from doy 1 00:30 and finish doy 365 or 366 at 23:30 


#PJG 2012  (Leap Year):
FLUX_SOIL<-read.table("C:/Users/Laura Morillas/Documents/data/PJ_guirdled/soil/From bandage script/PJ_girdle_2012_corrected_tstamps_23x_fluxall.dat", sep = "\t", header = TRUE) #sep = "\t" means separated by Tab
rawData<-FLUX_SOIL[,161:ncol(FLUX_SOIL)]   #to have only data coming from CR23X (FROM COLUMN 162 OF THE ORIGINAL FILE) including data doy 1 00:30 and finish doy 365 or 366 at 23:30 


fluxAll <- read.table("C:/Users/Laura Morillas/Documents/data/PJ_guirdled/Flux_all files/PJ_girdle_FLUX_all_2012.txt",sep = "\t", header = TRUE)
fluxAll[,ncol(fluxAll)]<-FLUX_SOIL$rain_Tot     #Because Rain_tot in  PJ_girdle_FLUX_all_2012.txt is wrong

PJC_rawData <- read.table("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/PJC_soil_complete/PJ_2012_soil_complete_07242014.txt", sep = ",", header = TRUE)
PJC_rawData<-PJC_rawData[2:17568,]            #only for 2011 data because:all soil_Year_complete datasets begin doy 1 00:30 and finish doy 365 or 366 at 23:30


                                        #all Flux_all files begin doy 1 00:00 and finish doy 365 or 366 at 23:30
                                        

### Select only the columns that contain SWC content. Searching here for *WC* as output
### by the CR23x. Also include the matlab datenum timestamp for POSIX conversion
rawData <- rawData[, c(1, grep("WC", names(rawData)))] #grep functions look for the colnames that includes "WC"
#PJC_rawData<-PJC_rawData[, c(1, grep("WC", names(PJC_rawData)))]

### Include precipitation in the working dataset (from flux all files): 
#rawData$PRECIP<-fluxAll[2:nrow(fluxAll),ncol(fluxAll)]   #This will work always because soil_Year_complete datasets begin doy 1 00:30 and finish doy 365 or 366 at 23:30 
                                                          #and Flux_all files begin doy 1 00:00 and finish doy 365 or 366 at 23:30. Also Precipitation is always the last column of FluxAll file
                                                          
                                                          
#Convert a numeric  MATLAB datenum (days since 0000-1-1 00:00) to seconds in 
#the Unix epoch (seconds since 1970-1-1 00:00). Specify a time zone if the 
#input datenum is anything other than the GMT/UTC time zone. 
# -- From the R file exchange --
matlab2POS = function(x, timez = "UTC") {
	days = x - 719529 	# 719529 = days from 1-1-0000 to 1-1-1970
	secs = days * 86400 # 86400 seconds in a day
	# This next string of functions is a complete disaster, but it works.
	# It tries to outsmart R by converting the secs value to a POSIXct value
	# in the UTC time zone, then converts that to a time/date string that 
	# should lose the time zone, and then it performs a second as.POSIXct()
	# conversion on the time/date string to get a POSIXct value in the user's 
	# specified timezone. Time zones are a goddamned nightmare.
	return(as.POSIXlt(strftime(as.POSIXct(secs, origin = '1970-1-1', 
			tz = 'UTC'), format = '%Y-%m-%d %H:%M', 
			tz = 'UTC', usetz = FALSE), tz = timez))
}


#For PJC 2011 :
timestamp <- matlab2POS(rawData$timestamps)
timestamp <- matlab2POS(PJC_rawData$timestamps)
#timestamp <- matlab2POS(rawData$tstamps)           # because timestamps in this datafile is  "tstamps"    matlab2POS doesn't work with PJG timestamps from the file "PJ_girdle_2011_corrected_tstamps_23x_fluxall.dat" at PJG 
                                                    #test_timestamp<-matlab2POS(Fluxsoil$timestamp)      #This is exaclty wrong than timestamp <- matlab2POS(rawData$tstamps) (using the 2 first columns of"PJ_girdle_2011_corrected_tstamps_23x_fluxall.dat")
                                                    #test2_timestamp<-matlab2POS(Fluxsoil$timestamp2)    #This is exaclty wrong than timestamp <- matlab2POS(rawData$tstamps) (using the 2 first columns of"PJ_girdle_2011_corrected_tstamps_23x_fluxall.dat")                                       
#For PJG 2011 ( UGLY BANDAGE):
PJC_timestamp <- matlab2POS(PJC_rawData$timestamps)
timestamp <- PJC_timestamp 
rawData$tstamps <- PJC_rawData$timestamps



timestampDF <- data.frame(ts = timestamp, year =  timestamp$year + 1900, month = timestamp$mon+1,
	mday = timestamp$mday, doy = (timestamp$yday)+1, 
	hrmin = round((timestamp$hour + timestamp$min / 60),digits = 1))
	

##########################################
# CREATING OUTPUT FILE
#
#Creating a matrix to record output files 
##########################################


OUTPUT0<- mat.or.vec(nrow(timestampDF), nc=length(names(rawData)))
OUTPUT<-data.frame(OUTPUT0)
names(OUTPUT)<-names(rawData[,c(1:dim(rawData)[2])])
OUTPUT[,1]<-rawData[,1]

GAP_account0<- mat.or.vec(nr=6, nc=length(names(rawData)))
GAP_account<-data.frame(GAP_account0)
names(GAP_account)<- c("Gap_level",names(rawData[,c(2:dim(rawData)[2])]))
GAP_account[,1]<-c("Raw without out of bounds","Gaps 1.5 hours filled","After Med Filt","After refining","Daily gaps from refined","Dry periods gap filled")





#######################################################
#######################################################
#
# FIRST STAGE: MEDIAN FILTER and general thresholds
#
#######################################################
#######################################################

########################
### CURRENT LOOP START
for(j in 1:ncol(rawData)){

colname <- names(rawData)[j+1]           #To determine the probe that we are going to check
currentCol <- data.frame(rawData[, j+1])
names(currentCol) <- eval(colname)

gapdata <- mat.or.vec(nrow(timestampDF), 1)   #mat.or.vec(nr,nc) creates and a zero vector of length nr if nc equals 1.This vectros will show the remaining gaps along the filtering process 


####################STEP 1 (REMOVING OUT OF BOUNDS VALUES)###########################################
## Check for bad values (not yet NaN) and correct
## Given these are soil water data, set all values
## greater than 1 and less than 0 to NaN            
## Need to establish a flag here for out of bounds values
currentCol[currentCol > 1] <- NaN
currentCol[currentCol < 0] <- NaN

## Count the gaps in the timeseries
for( i in 1:nrow(timestampDF)){
    if(is.na(currentCol[i,1])){
	gapdata[i] <- timestampDF$ts[i] 
    }else{
	gapdata[i] <- NaN 
    }
}

## Append the gap data to the dataframe
SWCdf <- data.frame(timestampDF,
	SWC = currentCol[,1] , GAP = gapdata)

### Reshape the data to long format for analysis / plotting
SWCdf_m <- melt(SWCdf, c('ts','year','month','mday','doy','hrmin','GAP'))

### Start the plotting by creating a simple theme to cut the frass                                   
theme_publish <- function (base_size = 12, base_family = "")                                      
{                                                                                                 
    theme_bw(base_size = base_size, base_family = base_family) %+replace%                         
        theme(legend.background = element_blank(), legend.key = element_blank(),                  
            panel.background = element_blank(),                                                   
            strip.background = element_blank(), plot.background = element_blank())                
} 

### Create a quick plot of the data for reference

ggplot(SWCdf_m, aes(ts, value)) + geom_line(aes(color = variable)) + 
	facet_grid(variable~., scales = 'free') + scale_color_manual(labels = c('GAPS','Soil Water'),
	values = c('red','black')) + geom_vline(aes(xintercept = GAP, alpha = .7, color = 'red'), SWCdf) +
	xlab('Time') + ylab('Measurements') + ggtitle('Raw SWC Timeseries') + theme_publish()
 
  GAP_n<-myfx(gapdata)
  
####################STEP 2 (LINEARLY GAP FILLING small, 1 hr or less, GAPS )###########################################
### If there is a leading gap thats small ( 1 hr or less) fill it with the 
### data that follows it. We can generalize this easily to handle larger
### leading gaps.
if(is.na(SWCdf$SWC[1]) & !is.na(SWCdf$SWC[2])){
SWCdf$SWC[1] <- SWCdf$SWC[2]
}

### Fill small gaps (1 hr or less) using na.approx and a linear function
filled_1 <- na.approx(SWCdf$SWC, SWCdf$ts, maxgap = 2, na.rm = FALSE)
gapdata_1 <- mat.or.vec(nrow(timestampDF), 1)
                                                                
                                                                
## Count the gaps in the timeseries
for( i in 1:nrow(timestampDF)){
    if(is.na(filled_1[i])){
	gapdata_1[i] <- timestampDF$ts[i] 
    }else{
	gapdata_1[i] <- NaN 
    }
}

## Append the gap data
SWCdf_1 <- data.frame(timestampDF, 
	SWC_fill1 = filled_1, GAP = gapdata_1)

## Melt the data to long format
SWCdf_1_m <- melt(SWCdf_1, c('ts','year','month','mday','doy','hrmin','GAP'))

ggplot(SWCdf_1_m, aes(ts, value)) + geom_line(aes(color = variable)) + 
	facet_grid(variable~., scales = 'free') + scale_color_manual(labels = c('GAPS','Soil Water'),
	values = c('red','black')) + geom_vline(aes(xintercept = GAP, alpha = .7, color = 'red'), SWCdf_1) +
	xlab('Time') + ylab('Measurements') + ggtitle('First Filtering SWC Timeseries') + theme_publish()
	
	 GAP_n_1<-myfx(gapdata_1)
	 
####################STEP 3 (APPLYING MEIAN FILTER) ###########################################
### Apply a median filter to the small gaps filled SWC series. We do this before gapfilling to avoid
### fitting curves to bad (not real) trends
### The window size (30) was determined to be the best in general cases, but if 
### the input data are difficult to smooth, changing this parameter (normally reducing the size)
### may help.

medfilteredSWC <- med.filter(SWCdf_1$SWC_fill1,50)
plot(medfilteredSWC)
filled_2<- unlist(medfilteredSWC$level)

gapdata_2 <- mat.or.vec(nrow(timestampDF), 1)

for( i in 1:nrow(timestampDF)){
    if(is.na(filled_2[i])){
	gapdata_2[i] <- timestampDF$ts[i] 
    }else{
	gapdata_2[i] <- NaN 
    }
}                                         #aT THIS POINT THIS gapdata_2 is recording the gaps in the previous filled dataset (SWCdf_1) instead the remaining gaps after the median filter 

      
        SWCdf_2 <- data.frame(timestampDF,
	SWC_fill2 = unlist(medfilteredSWC$level), GAP = gapdata_2)

## Melt the data to long format
SWCdf_2_m <- melt(SWCdf_2, c('ts','year','month','mday','doy','hrmin','GAP'))

#plot
ggplot(SWCdf_2_m, aes(ts, value)) + geom_line(aes(color = variable)) + 
	facet_grid(variable~., scales = 'free') + scale_color_manual(labels = c('GAPS','Soil Water'),
	values = c('red','black')) + geom_vline(aes(xintercept = GAP, alpha = .7, color = 'red'), SWCdf_2) +
	xlab('Time') + ylab('Measurements') + ggtitle('Median filter SWC Timeseries') + theme_publish()

   GAP_n_2<-myfx(gapdata_2) 

OUTPUT[,1+j]<-SWCdf_2[,7]   ##in progress I is not 1 anymore.....i changes during the loop
GAP_account[1,1+j]<-GAP_n
GAP_account[2,1+j]<-GAP_n_1
GAP_account[3,1+j]<-GAP_n_2


}

OUTPUT<-cbind(timestampDF[,2:ncol(timestampDF)],OUTPUT)


#Saving results:
#PJC
setwd('C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS')

#PJG
#setwd('C:/Users/Laura Morillas/Documents/data/PJ_guirdled/soil/SWC_FILTERED OUTPUTS')

write.csv(OUTPUT,paste(OUTPUT1_name,'.CSV',sep=''),row.names=TRUE)
write.csv(GAP_account,paste(GAP_name,'.CSV',sep=''),row.names=TRUE)
          #write.csv(GAP_account,paste(Title2,'.CSV',sep=''),row.names=TRUE)


##########################
#Ploting results STAGE 1:
##########################
OUTPUT<-OUTPUT[,6:ncol(OUTPUT)]

#AUTOMATED LOOP TO SAVE THE diagnostic plots PLOTS
 dev.off()
 dev.off()

PLOT1_name
#PJC:                        
pdf(file = "C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS/PlotSWC_PJC_2013_Median_filtered_30min.pdf")
#PJG:
#pdf(file = "C:/Users/Laura Morillas/Documents/data/PJ_guirdled/soil/SWC_FILTERED OUTPUTS/PlotSWC_PJG_2009_Median_filtered_30min.pdf")

par(mfrow=c(3,1),mar=c(2, 4, 2, 0.5),oma=c(2,0,2,4))

for( i in 2:ncol(rawData)) {
        
plot(timestampDF[,5],rawData[,i],type="l",lwd=2,ylim=c(0,0.5),ylab="SWC",xlab="DOY",main=names(rawData)[i])
lines(timestampDF[,5],OUTPUT[,i],type="l",col="red",lwd=2)
par(new=TRUE)                                                                           #This is to create a secondary axis in the same plot
plot(timestampDF[,5],fluxAll[,ncol(fluxAll)],type="h",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,max(fluxAll[,ncol(fluxAll)],na.rm=T)))        #This is to include precipitation
Yincrease<-(max(fluxAll[,ncol(fluxAll)],na.rm=T))/4  
axis( side=4)
mtext("Precip",side=4,line=3,cex=0.8)                                                                           #4 means that the rain axis is going to be writen in the right hand of the plot
legend("topleft",col=c("black","red","blue"),lty=1,legend=c("raw SWC","Filtered SWC","Precipitation"))

     }
     
dev.off()
dev.off()
dev.off()


      ###########################################
      #Determining best window size for median filter :
      ############################################

#      OUTPUT_WS30<-OUTPUT
#       OUTPUT_WS50<-OUTPUT
#       OUTPUT_WS70<-OUTPUT
#        OUTPUT_WS100<-OUTPUT
#          OUTPUT_WS150<-OUTPUT
#

        #To choose specific time periods for the plot:
#        RA<-c(1:nrow(timestampDF))
#        timestampDF_col<-cbind(timestampDF,RA)
#        TE<-range(timestampDF_col[timestampDF_col$doy>=200 & timestampDF_col$doy<=270,7])

#        Tsubset<-c(TE[1]:TE[2])

#        pdf(file = "C:/Users/Laura Morillas/Documents/R/SWC filtering/Diagnostic plots/PJC_2011_SWCfilter_WScompDetail30to70_2.pdf")
#        par(mfrow=c(3,1),mar=c(2, 4, 2, 0.5),oma=c(2,0,2,4))

#        for( i in 2:ncol(rawData)) {
        
#        plot(timestampDF_col[Tsubset,7],rawData[Tsubset,i],pch=8,type="p",lwd=2,ylim=c(0,0.5),ylab="SWC",main=names(rawData)[i]) #xlab="DOY",
#        lines(timestampDF_col[Tsubset,7],OUTPUT_WS30[Tsubset,i],pch=4,cex=0.3 ,type="p",col="pink",lwd=2)
#        lines(timestampDF_col[Tsubset,7],OUTPUT_WS50[Tsubset,i],pch=4,cex=0.3 ,type="p",col="green",lwd=2)
#       lines(timestampDF_col[Tsubset,7],OUTPUT_WS70[Tsubset,i],pch=4,cex=0.3 ,type="p",col="purple",lwd=2)
#        lines(timestampDF_col[Tsubset,7],OUTPUT_WS100[Tsubset,i],pch=4,cex=0.3 ,type="p",col="orange",lwd=2)
#        lines(timestampDF_col[Tsubset,7],OUTPUT_WS150[Tsubset,i],pch=4,cex=0.3 ,type="p",col="red",lwd=2)
#        par(new=TRUE)                                                                           #This is to create a secondary axis in the same plot
#        plot(timestampDF_col[Tsubset,7],fluxAll[Tsubset,ncol(fluxAll)],type="h",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,max(fluxAll[,ncol(fluxAll)],na.rm=T)))        #This is to include precipitation
#        Yincrease<-(max(fluxAll[,ncol(fluxAll)],na.rm=T))/4  
#        axis( side=4)
#        mtext("Precip",side=4,line=3,cex=0.8)                                                                           #4 means that the rain axis is going to be writen in the right hand of the plot
#        legend("topright",col=c("black","pink","green","purple","orange","red","blue"),lty=1,legend=c("raw SWC","Filtered SWC WS30","Filtered SWC WS50"
#        ,"Filtered SWC WS70","Filtered SWC WS100","Filtered SWC WS150","Precipitation"))

                     } 
     
#        dev.off()
#        dev.off()
#        dev.off()

#Note a value >150 for the window size parameter in the mediam filter reduces too much the higher values of SWC as a result of rain.
# window size of 30 or 50 seems to be the better option allowing the SWC keep realistic under fast SWC increases as a result to rain.wE DECIDED TO USE A Window size of 50.





#######################################################
#######################################################
#
# SECOND STAGE: Refining SWC 30min (using QC indicators)
#
#######################################################
#######################################################


######################################################
#####################################################
### Bulding a QC indicator  after filter SWC (ws=50)
#####################################################
#####################################################

#To generate the second filtered dataset:

OUTPUT_NEW<- mat.or.vec(nrow(OUTPUT), nc=length(names(OUTPUT)))
OUTPUT2<-data.frame(OUTPUT_NEW)
names(OUTPUT2)<-names(OUTPUT)
OUTPUT2[,1]<-OUTPUT[,1]    #PJC_rawData$timestamps IN CASE THAT OUTPUT[,1] is a constant value or a daily value
                           #For PJC2009 :OUTPUT2[,1]<-PJC_rawData$timestamps


#If you already filtered some of the probes

#PJC
#OUTPUT2<-read.csv("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS/SWC_PJC_2010_REFINED_30min.csv",header=T,sep = ",",dec=".")
#OUTPUT2<-OUTPUT2[,2:ncol(OUTPUT2)]
#PJG
OUTPUT2<-read.csv("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS/SWC_PJC_2013_REFINED_30min.csv",header=T,sep = ",",dec=".")
OUTPUT2<-OUTPUT2[,2:ncol(OUTPUT2)]

#This check has to be done probe by probe

#TRYING the diferential

probe=23   #change this one by one  (FROM 2 TO 28) 
df<-OUTPUT[,c(1,probe)]
df$PRECIP<-fluxAll[,ncol(fluxAll)]
head(df)


QC_dif <-rep(NA,nrow(df))
QC_2<-rep(0,nrow(df))
QC_3<-rep(0,nrow(df))
QC_4<-rep(0,nrow(df))
QC_5<-rep(0,nrow(df))
QC_6<-rep(0,nrow(df))
LP<-rep(NA,nrow(df))  #LP=LAST PRECIPITATION

###one

for (i in 2:nrow(df)){

    QC_dif[i]<-((df[i,2] - df[(i-1),2])*100 )/df[(i-1),2]

                     }
###two                     
                     
for (i in 25:nrow(df)){ 
   
    LP[i]<-sum(df[c((i-24):i),3],na.rm=T)              #RAIN IN THE LAST HALF A DAY
  
  ifelse(is.na(LP[i])==T | is.na(QC_dif[i])==T,QC_2[i]<-0,ifelse(LP[i]>0 & QC_dif[i]>0,QC_2[i]<-0, QC_2[i]<-QC_dif[i] ))  

#  if(is.na(LP[i])==T | is.na(QC_dif[i])==T) {QC_2[i]<-0}
#  if( LP[i]>0 & QC_dif[i]>0){   QC_2[i]<-0 }       #if it rained in the last half a day and SWC is increasing
    
#  else {QC_2[i]<-QC_dif[i]}
                      } 
                      
### three

#To determine the maximum real drying speed: 
df$row<-c(1:nrow(df))
df$QC_dif<-QC_dif

#To select the drying period we trust
dev.off()
Tsubset<-c(1:nrow(df))   #Tsubset has to be determine manually for the moment
plot(df[Tsubset,1],df[Tsubset,2],main=names(df)[2],ylim=c(0,0.5))
par(new=TRUE)
plot(df[Tsubset,1],fluxAll[Tsubset,ncol(fluxAll)],type="h",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,max(fluxAll[,ncol(fluxAll)],na.rm=T)))


#To detect the observation when maximum rain occurs:
T_Pmax<-max(df[df$PRECIP==max(df$PRECIP,na.rm=T),4],na.rm=T)  #To find when the maximun reduction of SWC after the maximum rain  happened                     
T_Pmax
Tsubset<-c(5200:14800)   #Tsubset has to be determine manually for the moment
plot(df[Tsubset,1],df[Tsubset,2],main=names(df)[2])
par(new=TRUE)
plot(df[Tsubset,1],fluxAll[Tsubset,ncol(fluxAll)],type="h",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,max(fluxAll[,ncol(fluxAll)],na.rm=T)))        #This is to include precipitation

Drymax<-min(df[Tsubset,5],na.rm=T)  # -0.5560834, to determine the maximun velocity of dry     
Drymax
                                  
for (i in 48:nrow(df)){ 
   
    LP[i]<-sum(df[c((i-47):i),3],na.rm=T)
    
  if(is.na(LP[i])==T | is.na(QC_dif[i])==T) {QC_3[i]<-0}
  if( LP[i]==0 & QC_2[i]>0){   QC_3[i]<-1 }                 # IF it didn't rain in the last day and SWC is increasing->BAD
  if(QC_2[i]< Drymax){   QC_3[i]<-1 }                      # the SWC can never be dry faster than the real maximum drying speed (Drymax), if it is->BAD

  else {QC_3[i]<-0}                                             #otherwise is OK
                      }
                      

####four

df$QC_3<-QC_3

windowsize = 32                           #half a day WINDOW

Q<-round(windowsize)                  #THE FIRST ROW INCLUDED IN WINDOW
Z<-(round(nrow(df)-windowsize/2)-1)   #THE LAST ROW INCLUDED IN WINDOW

#Depurating QC_indicator WITH STANDARD DEVIATION

for (i in Q:Z){

window <- df[(i - windowsize):i, ]
      
QC_4[i] <- sd(window[,6]) #use var() or standard error over QC_3

}

QC_5<-ifelse( QC_4 >=(0.49*max(QC_4,na.rm=T)),1,0)
QC_6<-ifelse( QC_4 >0,1,0)

df$QC_6<-QC_6
df$QC_5<-QC_5


#TO CHECK QC indicators


Tsubset<-c(1:nrow(df))  #For all data
Tsubset<-c(11000:13007)  #For monsoon
Tsubset<-c((T_Pmax-10):11500)  #For monsoon reduced
Tsubset<-c(5000:7000) #For spring
Tsubset<-c(1000:4500) #For beginning of year
Tsubset<-c(14000:17519) #For the year's end
Tsubset<-c(16150:16270)
Tsubset<-c(1000:4000) #For beginning of year
Tsubset<-c(7000:8000)

par(mfrow=c(4,1))

plot(df[Tsubset,1],df[Tsubset,2],main=names(df)[2],col="red",ylim=c(0,0.5))
par(new=TRUE)
plot(df[Tsubset,1],rawData[Tsubset,probe],type="l",lwd=2,ylim=c(0,0.5))
par(new=TRUE)
plot(df[Tsubset,1],QC_dif[Tsubset],ylim=c(0,30),col="orange",type="h")#,ylim=c(0,0.005)
par(new=TRUE)
plot(df[Tsubset,1],QC_3[Tsubset],ylim=c(0,1),col="purple",type="h")#,ylim=c(0,0.005)
par(new=TRUE)
plot(df[Tsubset,1],df[Tsubset,2],main=names(df)[2],,col="red",ylim=c(0,0.5))


plot(df[Tsubset,1],df[Tsubset,2],main=names(df)[2])
par(new=TRUE)
plot(df[Tsubset,1],QC_4[Tsubset],col="red",type="h")#,ylim=c(0,0.005)
par(new=TRUE)
plot(df[Tsubset,1],df[Tsubset,2],main=names(df)[2])

plot(df[Tsubset,1],df[Tsubset,2],main=names(df)[2])
par(new=TRUE)
plot(df[Tsubset,1],QC_5[Tsubset],col="brown",type="h")#,ylim=c(0,0.005)
par(new=TRUE)
plot(df[Tsubset,1],df[Tsubset,2],main=names(df)[2])
par(new=TRUE)
plot(df[Tsubset,1],fluxAll[Tsubset,ncol(fluxAll)],type="h",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,max(fluxAll[2:nrow(fluxAll),ncol(fluxAll)],na.rm=T)))        #This is to include precipitation


plot(df[Tsubset,1],df[Tsubset,2],main=names(df)[2])
par(new=TRUE)
plot(df[Tsubset,1],QC_6[Tsubset],col="pink",type="h")#,ylim=c(0,0.005)
par(new=TRUE)
plot(df[Tsubset,1],df[Tsubset,2],main=names(df)[2])

par(new=TRUE)
plot(df[Tsubset,1],fluxAll[Tsubset,ncol(fluxAll)],type="h",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,max(fluxAll[2:nrow(fluxAll),ncol(fluxAll)],na.rm=T)))        #This is to include precipitation

dev.off()

###################################################
#CLEANING THE SWC IN CONSIDERATION OF qc INDICATOR
###################################################

#Here depending on the plots we can chose using one of the
# Qc INDICATORs (VAR_DF_X) to filter or even none

####No extra cleaning:
        df$SWC_clean<-df[,2]
    
    
#####Cleaning:

    df$SWC_clean<-ifelse(df$QC_5==1,NA,df[,2])
    #df$SWC_clean2<-ifelse(df$QC_6==1,NA,df[,2])  #For comparison using diferent QC indicator

    #USING DIFERENT QC INDICATOR DEPENDING ON THE PERIOD
    df$SWC_clean<-ifelse(df$QC_5==1,NA,df[,2])
    df$SWC_clean6<-ifelse(df$QC_6==1,NA,df[,2])
    
    df$SWC_clean<-ifelse(df$row<T_Pmax ,df$SWC_clean6,df[,2])
    df$SWC_clean<-ifelse(df$row<4000,NA,df$SWC_clean)  #df$row<4000
    df$SWC_clean<-ifelse(df$row>3254 & df$row<3277,  NA,df$SWC_clean)
    
    df$SWC_clean<-ifelse(df$row>15000 ,df$SWC_clean6,df$SWC_clean)#df$SWC_clean
    df$SWC_clean<-ifelse(df$row>16158 & df$row<16226,df$SWC_clean5,df$SWC_clean)#df$SWC_clean
    
    df$SWC_clean<-ifelse(df$QC_6==1,NA,df[,2])
    df$SWC_clean<-ifelse(df$row>2169 & df$row<2249 & df$SWC_clean<=0.17,NA,df$SWC_clean)
    df$SWC_clean<-ifelse(df$row>16030 & df$row<16038,NA,df$SWC_clean)
    
    df$SWC_clean<-ifelse(df$row>15570 & df$row<15850 & df$SWC_clean<0.25,NA,df$SWC_clean)
    
#####BAD PROBE:
    df$SWC_clean<-rep(NA,nrow(df))
    
#####MANUAL CLEANING:
    min(df[df$row>T_Pmax & df$QC_6>0,4])
    
    df[df$row== 11000),1]
    abline(v= df[df$row==10949,1],col="green")
    df$SWC_clean<-ifelse(df$row>12230 & df$row<15375,NA,df[,2])  #11037
    
#######   When maximum daily values are ok but there is noise

#adding doy to df
df$doy<-timestampDF$doy

#Computing daily maximum SWC
SWC_MAX<-aggregate(df[,2], by=list(df$doy),FUN=max, na.rm=TRUE) 


#adding a column on df for SWC MAX
for (i in 1:nrow(df)){
df$SWC_MAX[i]<-SWC_MAX[SWC_MAX$Group.1==df$doy[i],2]
                     }
#cleaning one period                     
df$SWC_clean<-ifelse(df$QC_5==1 & df$SWC_clean<df$SWC_MAX ,NA,df$SWC_clean)
df$SWC_clean<-ifelse(df$row>13500 & df$row<13630 & df$SWC_clean<0.21 ,NA,df$SWC_clean)
     

#To check clean SWC

par(mfrow=c(3,1))

Tsubset<-c(1:nrow(df))  #For all data
Tsubset<-c(12000:14000)  #For monsoon
Tsubset<-c(5000:7000) #For spring
Tsubset<-c(1:5000) #For beginning of year
Tsubset<-c(15000:17519) #For the year's end

plot(df[Tsubset,1],df[Tsubset,2],ylim=c(0,0.5),main=names(df)[2])
par(new=TRUE)
#plot(df[Tsubset,1],df$QC_dif[Tsubset],,ylim=c(0,0.5),main=names(df)[2],col="pink")
#par(new=TRUE)
#plot(df[Tsubset,1],rawData[Tsubset,probe],ylim=c(0,0.5),main=names(df)[2],col='orange')
par(new=TRUE)   
plot(df[Tsubset,1],df$SWC_clean[Tsubset],ylim=c(0,0.5),main=names(df)[2],col="brown")
par(new=TRUE)

plot(df[Tsubset,1],fluxAll[Tsubset,ncol(fluxAll)],type="h",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,max(fluxAll[,ncol(fluxAll)],na.rm=T)))        #This is to include precipitation


#######################################
### Saving the clean SWC probe by probe
#OUTPUT2<-read.csv("C:/Users/Laura Morillas/Documents/data/PJ_guirdled/soil/SWC_FILTERED OUTPUTS/SWC_PJG_2010_REFINED_30min.csv",header=T,sep = ",",dec=".")

names(df)[2]

OUTPUT2$WC_O2_5_AVG<-df$SWC_clean     #Change here manually the name of the probe

head(OUTPUT2)


#################################
#Saving the filtered_Stage2 data:
#PJC:
setwd('C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS')
write.csv(OUTPUT2,paste(OUTPUT2_name,'.CSV',sep=''),row.names=TRUE)



#PJG:
#setwd('C:/Users/Laura Morillas/Documents/data/PJ_guirdled/soil/SWC_FILTERED OUTPUTS')
#write.csv(OUTPUT2,paste(OUTPUT2_name,'.CSV',sep=''),row.names=TRUE)

#aDDING TIME COLUMNS

OUTPUT2<-cbind(timestampDF,OUTPUT2[,2:ncol(OUTPUT2)])

#PJG:
setwd('C:/Users/Laura Morillas/Documents/data/PJ_guirdled/soil/SWC_FILTERED OUTPUTS')
write.csv(OUTPUT2,paste(OUTPUT2_name,'.CSV',sep=''),row.names=TRUE)

###################################################
#Accounting the new gaps (After refining30min SWC):  
  
#PJC:
#OUTPUT2<-read.csv("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS/SWC_PJC_2010_REFINED_30min.csv",header=T,sep = ",",dec=".")
#OUTPUT2<-OUTPUT2[,2:ncol(OUTPUT2)]

#PJG:
OUTPUT2<-read.csv("C:/Users/Laura Morillas/Documents/data/PJ_guirdled/soil/SWC_FILTERED OUTPUTS/SWC_PJG_2009_REFINED_30min.csv",header=T,sep = ",",dec=".")



OUTPUT2_sel<-OUTPUT2[,6:ncol(OUTPUT2)]


Gap_3 <- rep(NA, ncol(rawData[,2:ncol(rawData)]))

for( i in grep("WC",names(OUTPUT2_sel))){
     Gap_3[i-1]<- nrow(OUTPUT2_sel)-myfx(OUTPUT2_sel[,i])
    }

    
GAP_account[4,2:ncol(GAP_account)]<- Gap_3

GAP_account

#PJC:
setwd('C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS')
write.csv(GAP_account,paste(GAP_name,'.CSV',sep=''),row.names=TRUE)

#PJG:
#setwd('C:/Users/Laura Morillas/Documents/data/PJ_guirdled/soil/SWC_FILTERED OUTPUTS')
#write.csv(GAP_account,paste(GAP_name,'.CSV',sep=''),row.names=TRUE)
           
            #########################################################
            #charging data after second filtering stage if necessary:
            #########################################################

            rawData<- read.table("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/PJC_soil_complete/PJ_2011_soil_complete_07242014.txt", sep = ",", header = TRUE)
            rawData<-rawData[2:17520,]            #only for 2011 data because:all soil_Year_complete datasets begin doy 1 00:30 and finish doy 365 or 366 at 23:30
            rawData <- rawData[, c(1, grep("WC", names(rawData)))]

            OUTPUT<-read.csv("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS/SWC_PJC_2013_Median_filtered_30min.csv",header=T,sep = ",",dec=".")
            OUTPUT<-OUTPUT[,2:ncol(OUTPUT)]

            OUTPUT2<-read.csv("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS/SWC_PJC_2013_REFINED_30min.csv",header=T,sep = ",",dec=".")
            OUTPUT2<-OUTPUT2[,2:ncol(OUTPUT2)]
            OUTPUT2_sel<-OUTPUT2[,7:ncol(OUTPUT2)]

            #GAP_account<-read.csv("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS/SWC_GAPS_PJC_2011_30min.csv",header=T,sep = ",",dec=".")
            GAP_account<-read.csv("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS/GAPS_SWC_PJC_2013.csv",header=T,sep = ",",dec=".")
            GAP_account<-GAP_account[,2:ncol(GAP_account)]

            fluxAll<- read.csv("C:/Users/Laura Morillas/Documents/data/PJ_control/Flux_all files/PJ_FLUX_all_2010.csv",header=T,sep = ",",dec=".")

    

#############################################
#############################################
# PLOTING RESULTS FROM STAGE 2
#   (Refined) OUTPUT2
############################################
############################################

#AUTOMATED LOOP TO SAVE THE diagnostic plots PLOTS
 dev.off()
 
PLOT2_name 
 
#PJC                      
pdf(file = "C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS/PlotSWC_PJC_2013_REFINED_30min.pdf")

#PJG                      
#pdf(file = "C:/Users/Laura Morillas/Documents/data/PJ_guirdled/soil/SWC_FILTERED OUTPUTS/PlotSWC_PJG_2012_REFINED_30min.pdf")
par(mfrow=c(3,1),mar=c(2, 4, 2, 0.5),oma=c(2,0,2,4))

for( i in 2:ncol(rawData)) {
        
plot(timestampDF[,5],rawData[,i],type="l",lwd=2,ylim=c(0,0.5),ylab="SWC",xlab="DOY",main=names(rawData)[i])
lines(timestampDF[,5],OUTPUT[,i],type="l",col="red",lwd=2)
lines(timestampDF[,5],OUTPUT2_sel[,i],type="l",col="green",lwd=2)
par(new=TRUE)                                                                           #This is to create a secondary axis in the same plot
plot(timestampDF[,5],fluxAll[,ncol(fluxAll)],type="h",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,max(fluxAll[,ncol(fluxAll)],na.rm=T)))        #This is to include precipitation
Yincrease<-(max(fluxAll[,ncol(fluxAll)],na.rm=T))/4  
axis( side=4)
mtext("Precip",side=4,line=3,cex=0.8)                                                                           #4 means that the rain axis is going to be writen in the right hand of the plot
legend("topleft",col=c("black","red","green","blue"),lty=1,legend=c("raw SWC","First_Filtered SWC","Second_Filtered SWC","Precipitation"))

     }

dev.off()
dev.off()
dev.off()





#######################################################
#######################################################
#
# THIRD STAGE: Estimating DAILY values and gap filling
#
#######################################################
#######################################################

###Estimating Daily values:

OUTPUT1_daily_avg<- aggregate(OUTPUT2, by=list(OUTPUT2$doy),FUN=mean, na.rm=TRUE)

fluxAll$doy<-as.integer(fluxAll$jday)
Daily_rain<-aggregate(fluxAll[,(ncol(fluxAll)-1):ncol(fluxAll)], by=list(fluxAll[,ncol(fluxAll)]),FUN=sum, na.rm=TRUE)

OUTPUT1_daily_avg$rain_Tot<-Daily_rain$rain_Tot
OUTPUT1_daily_avg<-OUTPUT1_daily_avg[,2:ncol(OUTPUT1_daily_avg)]

#Saving Daily SWC values:
#PJC:
setwd('C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS')
write.csv(OUTPUT1_daily_avg,paste(OUTPUT3_name,'.CSV',sep=''),row.names=TRUE)

#PJG:
#setwd('C:/Users/Laura Morillas/Documents/data/PJ_guirdled/soil/SWC_FILTERED OUTPUTS')
#write.csv(OUTPUT1_daily_avg,paste(OUTPUT3_name,'.CSV',sep=''),row.names=TRUE)

#########################
#Accounting the new gaps (daily scale):    
OUTPUT1_daily_avg_sel<-OUTPUT1_daily_avg[,6:ncol(OUTPUT1_daily_avg)]

Gap_4 <- rep(NA, ncol(rawData[,2:ncol(rawData)]))

for( i in grep("WC",names(OUTPUT1_daily_avg_sel))){
     Gap_4[i-1]<- nrow(OUTPUT1_daily_avg_sel)-myfx(OUTPUT1_daily_avg_sel[,i])
    }

    
GAP_account[5,2:ncol(GAP_account)]<- Gap_4



#PJC:
setwd('C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS')
write.csv(GAP_account,paste(GAP_name,'.CSV',sep=''),row.names=TRUE)

#PJG:
#setwd('C:/Users/Laura Morillas/Documents/data/PJ_guirdled/soil/SWC_FILTERED OUTPUTS')
#write.csv(GAP_account,paste(GAP_name,'.CSV',sep=''),row.names=TRUE)

########################################
#      Ploting Result 1 THIRD STAGE 
#     (daily averages of refined SWC)
########################################
dev.off()

PLOT3_name
#PJC:
pdf(file = "C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS/PlotSWC_PJC_2013_REFINED_Daily.pdf")

#PJG:
#pdf(file = "C:/Users/Laura Morillas/Documents/data/PJ_guirdled/soil/SWC_FILTERED OUTPUTS/PlotSWC_PJG_2012_REFINED_Daily.pdf")

par(mfrow=c(3,1),mar=c(2, 4, 2, 0.5),oma=c(2,0,2,4))

for( i in grep("WC",names(OUTPUT1_daily_avg))) {
        
plot(OUTPUT1_daily_avg$doy,OUTPUT1_daily_avg[,i],type="l",lwd=2,ylim=c(0,0.5),ylab="SWC",xlab="DOY",main=names(OUTPUT1_daily_avg)[i],col="green")
par(new=TRUE)                                                                           #This is to create a secondary axis in the same plot
plot(OUTPUT1_daily_avg$doy,Daily_rain$rain_Tot ,type="h",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,max(Daily_rain$rain_Tot)))        #This is to include precipitation
Yincrease<-max(Daily_rain)/4  
axis( side=4)
mtext("Precip",side=4,line=3,cex=0.8)                                                                           #4 means that the rain axis is going to be writen in the right hand of the plot
legend("topleft",col=c("black","blue"),lty=1,legend=c("Second_Filtered daily SWC","Precipitation"))

     }

dev.off()
dev.off()
dev.off()



######################################################
######################################################
### Gap filling daily data when there was not rain:
######################################################
######################################################


INPUT<-OUTPUT1_daily_avg[,c(5,7:ncol(OUTPUT1_daily_avg))]     #A DATASET only including doy and SWC measurements and rain
head(INPUT)
   

Gap_start<-rep(NA,nrow(INPUT))
Gap_end<-rep(NA,nrow(INPUT))
Gap_start2<-rep(NA,nrow(INPUT))
Gap_end2<-rep(NA,nrow(INPUT))
Gap_Rain<-rep(NA,nrow(INPUT))
SWC_FILLED<-rep(NA,nrow(INPUT))
SWC_Gapfilled<- mat.or.vec(nr=nrow(INPUT),nc=1+length(grep("WC",names(INPUT))))
SWC_Gapfilled[,1]<-INPUT[,1]
SWC_Gapfilled<-data.frame(SWC_Gapfilled)
names(SWC_Gapfilled)<-names(INPUT[,1:ncol(INPUT)-1])


for( j in grep("WC",names(INPUT))){

Gap_id<-ifelse( is.na(INPUT[,j])==TRUE,1,0)


Gap_start[1]<-ifelse(Gap_id[1]==1,INPUT[1,1],NA)
Gap_end[length(Gap_id)]<-ifelse(Gap_id[length(Gap_id)]==1,INPUT[length(Gap_id),1],NA)
                              

for (i in 2:(nrow(INPUT)-1)){

Gap_start[i]<-ifelse(Gap_id[i-1]==0 && Gap_id[i]==1,INPUT[i-1,1],NA)
Gap_end[i]<- ifelse(Gap_id[i+1]==0 && Gap_id[i]==1,INPUT[i+1,1],NA)
                              }
                              
for (i in 1:nrow(INPUT)){
                              
Gap_start2[i]<-ifelse(Gap_id[i]==1,max(Gap_start[1:i],na.rm=T),0)
Gap_end2[i]<- ifelse(Gap_id[i]==1,min(Gap_end[i:length(Gap_end)],na.rm=T),0)
Gap_Rain[i]<-ifelse(Gap_id[i]==1,sum(INPUT[Gap_start2[i]:Gap_end2[i],grep("rain_Tot",names(INPUT))],na.rm=T),NA)
                        }
                        
Gap_pos<- (Gap_start[is.na(Gap_start)!=T]+1)
                       
for (i in Gap_pos){

if(Gap_Rain[i]>0)  {SWC_FILLED[Gap_start2[i]:Gap_end2[i]]<-INPUT[Gap_start2[i]:Gap_end2[i],j]}
else {SWC_FILLED[Gap_start2[i]:Gap_end2[i]]<-na.approx(INPUT[Gap_start2[i]:Gap_end2[i],j])}
                        
                       }
                       
SWC_Gapfilled[,j]<-ifelse(Gap_id==1,SWC_FILLED,INPUT[,j])                       
                       
                          }
                          
SWC_Gapfilled$rain_Tot<-INPUT$rain_Tot

#PJC:
setwd('C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS')
write.csv(SWC_Gapfilled,paste(OUTPUT4_name,'.CSV',sep=''),row.names=TRUE)

#PJG:
setwd('C:/Users/Laura Morillas/Documents/data/PJ_guirdled/soil/SWC_FILTERED OUTPUTS')
write.csv(SWC_Gapfilled,paste(OUTPUT4_name,'.CSV',sep=''),row.names=TRUE)

########################################
#      Ploting Result 2 THIRD STAGE 
#   (daily refined and Gap filled SWC)
########################################

#AUTOMATED LOOP TO SAVE THE diagnostic plots PLOTS
 dev.off()
 
PLOT4_name
#PJC:
pdf(file = "C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS/PlotSWC_PJC_2013_Refined&Gapfilled_Daily.pdf")     #STICK HERE THE NAME OF plot 3

#PJG:
#pdf(file = "C:/Users/Laura Morillas/Documents/data/PJ_guirdled/soil/SWC_FILTERED OUTPUTS/PlotSWC_PJG_2012_Refined&Gapfilled_Daily.pdf")


par(mfrow=c(3,1),mar=c(2, 4, 2, 0.5),oma=c(2,0,2,4))

for( i in 2:ncol(rawData)) {
        
plot(INPUT$doy,SWC_Gapfilled[,i],type="l",lwd=2,ylim=c(0,0.5),ylab="SWC",xlab="DOY",main=names(rawData)[i],col="red")
lines(INPUT$doy,INPUT[,i],type="l",lwd=2)
par(new=TRUE)                                                                           #This is to create a secondary axis in the same plot
plot(INPUT$doy,INPUT$rain_Tot,type="h",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,max(INPUT$rain_Tot,na.rm=T)))        #This is to include precipitation
Yincrease<-(max(INPUT$rain_Tot,na.rm=T))/4  
axis( side=4)
mtext("Precip",side=4,line=3,cex=0.8)                                                                           #4 means that the rain axis is going to be writen in the right hand of the plot
legend("topleft",col=c("red","black","blue"),lty=1,legend=c("Gap filled SWC","Refined SWC","Precipitation"))

     }

dev.off()
dev.off()
dev.off()


#########################
#Accounting the new gaps (gap filled daily data):    
#Input:SWC_Gapfilled

Gap_5<- rep(NA, ncol(rawData[,2:ncol(rawData)]))

for( i in grep("WC",names(SWC_Gapfilled))){
     Gap_5[i-1]<- nrow(SWC_Gapfilled)-myfx(SWC_Gapfilled[,i])
    }

   
GAP_account[6,2:ncol(GAP_account)]<- Gap_5

    
#PJC:
setwd('C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS')
write.csv(GAP_account,paste(GAP_name,'.CSV',sep=''),row.names=TRUE)

#PJG:
setwd('C:/Users/Laura Morillas/Documents/data/PJ_guirdled/soil/SWC_FILTERED OUTPUTS')
write.csv(GAP_account,paste(GAP_name,'.CSV',sep=''),row.names=TRUE)








##############################################################################################
#
#INCLUDING AVERAGES BY LOCATION AND DEPTH (DAILY GAP FILLED)
#
#############################################################################################
#PJC:

SWC_PJC<-read.csv("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS/SWC_PJC_2012_Refined&Gapfilled_Daily.csv",header=T,sep = ",",dec=".")
names(SWC_PJC)

P_5<-c(3,6,9)
P_10<-c(4,7,10)
P_30<-c(5,8,11)

J_5<-c(12,15,18)
J_10<-c(13,16,19)
J_30<-c(14,17,20)

O_5<-c(21,24,27)
O_10<-c(22,25,28)
O_30<-c(23,26,29)

SWC_PJC$P_5<-rowMeans(SWC_PJC[,P_5], na.rm = T)
SWC_PJC$P_10<-rowMeans(SWC_PJC[,P_10], na.rm = T)
SWC_PJC$P_30<-rowMeans(SWC_PJC[,P_30], na.rm = T)

SWC_PJC$J_5<-rowMeans(SWC_PJC[,J_5], na.rm = T)
SWC_PJC$J_10<-rowMeans(SWC_PJC[,J_10], na.rm = T)
SWC_PJC$J_30<-rowMeans(SWC_PJC[,J_30], na.rm = T)

SWC_PJC$O_5<-rowMeans(SWC_PJC[,O_5], na.rm = T)
SWC_PJC$O_10<-rowMeans(SWC_PJC[,O_10], na.rm = T)
SWC_PJC$O_30<-rowMeans(SWC_PJC[,O_30], na.rm = T)

SWC_PJC$P_AVG<-rowMeans(SWC_PJC[,grep("P_",names(SWC_PJC))], na.rm = F)
SWC_PJC$J_AVG<-rowMeans(SWC_PJC[,grep("J_",names(SWC_PJC))], na.rm = F)
SWC_PJC$O_AVG<-rowMeans(SWC_PJC[,grep("O_",names(SWC_PJC))], na.rm = F)

SWC_PJC$P_AVG<-rowMeans(SWC_PJC[,grep("P_",names(SWC_PJC))], na.rm = F)
SWC_PJC$J_AVG<-rowMeans(SWC_PJC[,grep("J_",names(SWC_PJC))], na.rm = F)
SWC_PJC$O_AVG<-rowMeans(SWC_PJC[,grep("O_",names(SWC_PJC))], na.rm = F)

setwd('C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS')
write.csv(SWC_PJC,paste(OUTPUT5_name,'.CSV',sep=''),row.names=TRUE)


#PJG
SWC_PJG<-read.csv("C:/Users/Laura Morillas/Documents/data/PJ_guirdled/soil/SWC_FILTERED OUTPUTS/SWC_PJG_2010_Refined&Gapfilled_Daily.csv",sep=",",header=TRUE,dec=".")
names(SWC_PJG)

#THERE ARE ISSUES IN THE OUTPUT NAMES AT THIS SITE BECAUSE SENSORS ARE NOT ACTUALLY LOCATED AT 5,10 AND 30,as at PJC
#THEY ARE LOCATED AT 5,10 AND AS MUCH DEEP AS YOU CAN GO ABOVE THE CALICHE LAYER, SO tHIS PART OF THE CODE 
#TRY TO ACCOUNT OR FIX THOSE CONFUSSIONS AND TRY TO GENERATE DATA COMPARABLE TO pjc SITE. 

P_5<-c(3,6,9)
P_10<-c(4,7,10)
P_20<-c(5,8)
P_30<-c(11) #P3_30 CM IS ACTUALLY LOCATED AT 32.5 CM (ALMOST 30CM)


J_5<-c(12,15,18)
J_10<-c(13,16,19)
J_20<-c(14,17)
J_30<-c(20)    #J3_30 CM IS ACTUALLY LOCATED AT 30-35 CM

O_5<-c(21,24,27)
O_10<-c(22,25,28)
O_20<-c(23,26)
O_30<-c(29)   

SWC_PJG$P_5<-rowMeans(SWC_PJG[,P_5], na.rm = T)
SWC_PJG$P_10<-rowMeans(SWC_PJG[,P_10], na.rm = T)
SWC_PJG$P_20<-rowMeans(SWC_PJG[,P_20], na.rm = T)
SWC_PJG$P_30<-SWC_PJG[,P_30]


SWC_PJG$J_5<-rowMeans(SWC_PJG[,J_5], na.rm = T)
SWC_PJG$J_10<-rowMeans(SWC_PJG[,J_10], na.rm = T)
SWC_PJG$J_20<-rowMeans(SWC_PJG[,J_20], na.rm = T)
SWC_PJG$J_30<-SWC_PJG[,J_30]

SWC_PJG$O_5<-rowMeans(SWC_PJG[,O_5], na.rm = T)
SWC_PJG$O_10<-rowMeans(SWC_PJG[,O_10], na.rm = T)
SWC_PJG$O_20<-rowMeans(SWC_PJG[,O_20], na.rm = T)
SWC_PJG$O_30<-SWC_PJG[,O_30]

SWC_PJG$P_AVG<-rowMeans(SWC_PJG[, grep("P_",names(SWC_PJG))], na.rm = F)
SWC_PJG$J_AVG<-rowMeans(SWC_PJG[,grep("J_",names(SWC_PJG))], na.rm = F)
SWC_PJG$O_AVG<-rowMeans(SWC_PJG[,grep("O_",names(SWC_PJG))], na.rm = F)

SWC_PJG$P_AVG<-rowMeans(SWC_PJG[,grep("P_",names(SWC_PJG))], na.rm = F)
SWC_PJG$J_AVG<-rowMeans(SWC_PJG[,grep("J_",names(SWC_PJG))], na.rm = F)
SWC_PJG$O_AVG<-rowMeans(SWC_PJG[,grep("O_",names(SWC_PJG))], na.rm = F)

setwd('C:/Users/Laura Morillas/Documents/data/PJ_guirdled/soil/SWC_FILTERED OUTPUTS')
write.csv(SWC_PJG,paste(OUTPUT5_name,'.CSV',sep=''),row.names=TRUE)

