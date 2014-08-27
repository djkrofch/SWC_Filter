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
3
#
#
# For inquiries about this script's use or function, contact the authors at:
#	Dan Krofcheck      krofcheck@gmail.com
# Laura Morillas     laura.morillas.gonzalez@gmail.com
#
#################################################################################
#################################################################################

#For naming OUTPUT files:

SITES<-c("PJC","PJG")

SITE.CODE<-1              #WRITE HERE THE CODE (1 for PJC, 2 for PJG) OF THE SITE TO PROCESS
YEAR<-"2011"              #WRITE HERE THE YEAR OF YEARS OF DATA TO PROCESS

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
                        
rawData <- read.table("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/PJC_soil_complete/PJ_2011_soil_complete_07242014.txt", sep = ",", header = TRUE)
#rawData <- read.table("C:/Users/David/Documents/Laura/SWC filtering/PJ_2011_soil_complete_07242014.txt", sep = ",", header = TRUE)

rawData<-rawData[2:17520,]            #only for 2011 data because:all soil_Year_complete datasets begin doy 1 00:30 and finish doy 365 or 366 at 23:30
fluxAll  <- read.csv("C:/Users/Laura Morillas/Documents/data/PJ_control/Flux_all files/PJ_FLUX_all_2011_Dan.csv",header=T,sep = ",",dec=".")
#fluxAll  <- read.csv("C:/Users/David/Documents/Laura/SWC filtering/PJ_FLUX_all_2011_Dan.csv",header=T,sep = ",",dec=".")

                                        #all Flux_all files begin doy 1 00:00 and finish doy 365 or 366 at 23:30
                                        
### Select only the columns that contain SWC content. Searching here for *WC* as output
### by the CR23x. Also include the matlab datenum timestamp for POSIX conversion
rawData <- rawData[, c(1, grep("WC", names(rawData)))] #grep functions look for the colnames that includes "WC"

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

timestamp <- matlab2POS(rawData$timestamps)

timestampDF <- data.frame(ts = timestamp, year =  timestamp$year + 1900, month = timestamp$mon+1,
	mday = timestamp$mday, doy = (timestamp$yday)+1, 
	hrmin = round((timestamp$hour + timestamp$min / 60),digits = 1))
	


#Creating a matrix to record output files 

OUTPUT0<- mat.or.vec(nrow(timestampDF), nc=length(names(rawData)))
OUTPUT<-data.frame(OUTPUT0)
names(OUTPUT)<-names(rawData[,c(1:dim(rawData)[2])])
OUTPUT[,1]<-rawData$timestamp

GAP_account0<- mat.or.vec(nr=6, nc=length(names(rawData)))
GAP_account<-data.frame(GAP_account0)
names(GAP_account)<- c("Gap_level",names(rawData[,c(2:dim(rawData)[2])]))
GAP_account[,1]<-c("Raw without out of bounds","Gaps 1.5 hours filled","After Med Filt","After refining","Daily gaps from refined","Dry periods gap filled")





#######################################################
#######################################################
#
# FIRDT STAGE: MEDIAN FILTER and general thresholds
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


#Saving results:
setwd('C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS')

write.csv(OUTPUT,paste(OUTPUT1_name,'.CSV',sep=''),row.names=TRUE)
write.csv(GAP_account,paste(GAP_name,'.CSV',sep=''),row.names=TRUE)
          #write.csv(GAP_account,paste(Title2,'.CSV',sep=''),row.names=TRUE)


##########################
#Ploting results STAGE 1:
##########################


#AUTOMATED LOOP TO SAVE THE diagnostic plots PLOTS
 dev.off()
 dev.off()
                        
pdf(file = "C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS/PlotSWC_PJC_2011_Median_filtered_30min.pdf")
par(mfrow=c(3,1),mar=c(2, 4, 2, 0.5),oma=c(2,0,2,4))

for( i in 2:ncol(rawData)) {
        
plot(timestampDF[,5],rawData[,i],type="l",lwd=2,ylim=c(0,0.5),ylab="SWC",xlab="DOY",main=names(rawData)[i])
lines(timestampDF[,5],OUTPUT[,i],type="l",col="red",lwd=2)
par(new=TRUE)                                                                           #This is to create a secondary axis in the same plot
plot(timestampDF[,5],fluxAll[2:nrow(fluxAll),ncol(fluxAll)],type="h",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,max(fluxAll[2:nrow(fluxAll),ncol(fluxAll)],na.rm=T)))        #This is to include precipitation
Yincrease<-(max(fluxAll[2:nrow(fluxAll),ncol(fluxAll)],na.rm=T))/4  
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
#        plot(timestampDF_col[Tsubset,7],fluxAll[Tsubset,ncol(fluxAll)],type="h",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,max(fluxAll[2:nrow(fluxAll),ncol(fluxAll)],na.rm=T)))        #This is to include precipitation
#        Yincrease<-(max(fluxAll[2:nrow(fluxAll),ncol(fluxAll)],na.rm=T))/4  
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
OUTPUT2[,1]<-OUTPUT[,1]

#If you already filtered some of the probes
OUTPUT2<-read.csv("C:/Users/Laura Morillas/Documents/R/SWC filtering/PJC_11_Filtered2SWC_WS50,.csv",header=T,sep = ",",dec=".")


#This check has to be done probe by probe

#TRYING the diferential

probe=7  #change this one by one  (FROM 2 TO 28) 
df<-OUTPUT[,c(1,probe)]
df$PRECIP<-fluxAll[2:nrow(fluxAll),ncol(fluxAll)]

var_df<-rep(NA,nrow(df))
var_df_2<-rep(0,nrow(df))
var_df_3<-rep(0,nrow(df))
var_df_4<-rep(0,nrow(df))
var_df_5<-rep(0,nrow(df))
var_df_6<-rep(0,nrow(df))
LP<-rep(NA,nrow(df))  #LP=LAST PRECIPITATION

###one

for (i in 2:nrow(df)){

    var_df[i]<-((df[i,2] - df[(i-1),2])*100 )/df[(i-1),2]

                     }
###two                     
                     
for (i in 25:nrow(df)){ 
   
    LP[i]<-sum(df[c((i-24):i),3],na.rm=T)              #RAIN IN THE LAST HALF A DAY
    
  if( LP[i]>0 & var_df[i]>0){   var_df_2[i]<-0 }       #if it rained in the last half day and SWC is increasing
    
  else {var_df_2[i]<-var_df[i]}
                      } 
                      
### three

#To determine the maximum real drying speed: 
df$row<-c(1:nrow(df))
df$var_df<-var_df

#To select the drying period we trust
Tsubset<-c(1:nrow(df))   #Tsubset has to be determine manually for the moment
plot(df[Tsubset,1],df[Tsubset,2],main=names(df)[2],ylim=c(0,0.5))
par(new=TRUE)
plot(df[Tsubset,1],fluxAll[Tsubset,ncol(fluxAll)],type="h",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,max(fluxAll[2:nrow(fluxAll),ncol(fluxAll)],na.rm=T)))        #This is to include precipitation

#To detect the observation when maximum rain occurs:
max(df[df$PRECIP==max(df$PRECIP,na.rm=T),4],na.rm=T)  #To find when the maximun reduction of SWC after the maximum rain  happened                     

Tsubset<-c(10299:10400)   #Tsubset has to be determine manually for the moment
plot(df[Tsubset,1],df[Tsubset,2],main=names(df)[2])
par(new=TRUE)
plot(df[Tsubset,1],fluxAll[Tsubset,ncol(fluxAll)],type="h",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,max(fluxAll[2:nrow(fluxAll),ncol(fluxAll)],na.rm=T)))        #This is to include precipitation

Drymax<-min(df[Tsubset,5],na.rm=T)  # -0.5560834, to determine the maximun velocity of dry     
Drymax
                                  
for (i in 48:nrow(df)){ 
   
    LP[i]<-sum(df[c((i-47):i),3],na.rm=T)
#  if( LP[i]>0 & var_df_2[i]< (-0.56)){   var_df_3[i]<-1 }          # IF it rained in the last day and SWC is decreasing faster than the maximum real decrease after the maximum rain->BAD
  if( LP[i]==0 & var_df_2[i]>0){   var_df_3[i]<-1 }                 # IF it didn't rain in the last day and SWC is increasing->BAD
  if(var_df_2[i]< Drymax){   var_df_3[i]<-1 }                      # the SWC can never be larger than the real maximum, if it is->BAD

  else {var_df_3[i]<-0}                                             #otherwise is OK
                      }
                      

####four

df$var_df_3<-var_df_3

windowsize = 32                           #half a day WINDOW

Q<-round(windowsize)                  #THE FIRST ROW INCLUDED IN WINDOW
Z<-(round(nrow(df)-windowsize/2)-1)   #THE LAST ROW INCLUDED IN WINDOW

#Depurating QC_indicator WITH STANDARD DEVIATION

for (i in Q:Z){

window <- df[(i - windowsize):i, ]
      
var_df_4[i] <- sd(window[,6]) #use var() or standard error over var_df_3

}

var_df_5<-ifelse( var_df_4 >=(0.49*max(var_df_4,na.rm=T)),1,0)
var_df_6<-ifelse( var_df_4 >0,1,0)

df$var_df_6<-var_df_6
df$var_df_5<-var_df_5


#TO CHECK QC indicators

Tsubset<-c(1:nrow(df))  #For all data
Tsubset<-c(9600:13007)  #For monsoon
Tsubset<-c(5000:7000) #For spring
Tsubset<-c(1:5000) #For beginning of year
Tsubset<-c(15000:17519) #For the year's end
Tsubset<-c(16848:17519)

par(mfrow=c(4,1))

plot(df[Tsubset,1],df[Tsubset,2],main=names(df)[2],col="red",ylim=c(0,0.5))
par(new=TRUE)
plot(df[Tsubset,1],rawData[Tsubset,probe],type="l",lwd=2,ylim=c(0,0.5))
par(new=TRUE)
plot(df[Tsubset,1],var_df[Tsubset],ylim=c(0,30),col="orange",type="h")#,ylim=c(0,0.005)
par(new=TRUE)
plot(df[Tsubset,1],var_df_3[Tsubset],ylim=c(0,1),col="purple",type="h")#,ylim=c(0,0.005)
par(new=TRUE)
plot(df[Tsubset,1],df[Tsubset,2],main=names(df)[2],,col="red",ylim=c(0,0.5))


plot(df[Tsubset,1],df[Tsubset,2],main=names(df)[2])
par(new=TRUE)
plot(df[Tsubset,1],var_df_4[Tsubset],col="red",type="h")#,ylim=c(0,0.005)
par(new=TRUE)
plot(df[Tsubset,1],df[Tsubset,2],main=names(df)[2])

plot(df[Tsubset,1],df[Tsubset,2],main=names(df)[2])
par(new=TRUE)
plot(df[Tsubset,1],var_df_5[Tsubset],col="brown",type="h")#,ylim=c(0,0.005)
par(new=TRUE)
plot(df[Tsubset,1],df[Tsubset,2],main=names(df)[2])
par(new=TRUE)
plot(df[Tsubset,1],fluxAll[Tsubset,ncol(fluxAll)],type="h",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,max(fluxAll[2:nrow(fluxAll),ncol(fluxAll)],na.rm=T)))        #This is to include precipitation


plot(df[Tsubset,1],df[Tsubset,2],main=names(df)[2])
par(new=TRUE)
plot(df[Tsubset,1],var_df_6[Tsubset],col="pink",type="h")#,ylim=c(0,0.005)
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

    df$SWC_clean<-ifelse(df$var_df_5==1,NA,df[,2])
    #df$SWC_clean2<-ifelse(df$var_df_6==1,NA,df[,2])  #For comparison using diferent QC indicator

    #USING DIFERENT QC INDICATOR DEPENDING ON THE PERIOD
    df$SWC_clean5<-ifelse(df$var_df_5==1,NA,df[,2])
    df$SWC_clean6<-ifelse(df$var_df_6==1,NA,df[,2])
    df$SWC_clean<-ifelse(df$row>10299,df$SWC_clean5,df$SWC_clean6)
  

#####BAD PROBE:
    df$SWC_clean<-rep(NA,nrow(df))

#To check clean SWC

par(mfrow=c(3,1))

Tsubset<-c(1:nrow(df))  #For all data
Tsubset<-c(9600:13007)  #For monsoon
Tsubset<-c(5000:7000) #For spring
Tsubset<-c(1:5000) #For beginning of year
Tsubset<-c(15000:17519) #For the year's end

plot(df[Tsubset,1],df[Tsubset,2],ylim=c(0,0.5),main=names(df)[2])
par(new=TRUE)
plot(df[Tsubset,1],df$SWC_clean[Tsubset],ylim=c(0,0.5),main=names(df)[2],col="brown")
#par(new=TRUE)
#plot(df[Tsubset,1],df$SWC_clean[Tsubset],,ylim=c(0,0.5),main=names(df)[2],col="pink")
par(new=TRUE)
plot(df[Tsubset,1],fluxAll[Tsubset,ncol(fluxAll)],type="h",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,max(fluxAll[2:nrow(fluxAll),ncol(fluxAll)],na.rm=T)))        #This is to include precipitation


#######################################
### Saving the clean SWC probe by probe

names(df)[2]

OUTPUT2$x37WC_P2_30_AVGH<-df$SWC_clean     #Change here manually the name of the probe

head(OUTPUT2)


#################################
#Saving the filtered_Stage2 data:
setwd('C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS')

write.csv(OUTPUT2,paste(OUTPUT2_name,'.CSV',sep=''),row.names=TRUE)
#write.csv(OUTPUT2,paste('PJC_11_Filtered2SWC_WS50','.CSV',sep=','),row.names=TRUE)



###################################################
#Accounting the new gaps (After refining30min SWC):  
  
OUTPUT2<-read.csv("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS/SWC_PJC_2011_REFINED_30min.csv",header=T,sep = ",",dec=".")
OUTPUT2<-OUTPUT2[,2:ncol(OUTPUT2)]
OUTPUT2_sel<-OUTPUT2[,6:ncol(OUTPUT2)]


Gap_3 <- rep(NA, ncol(rawData[,2:ncol(rawData)]))

for( i in grep("WC",names(OUTPUT2_sel))){
     Gap_3[i-1]<- nrow(OUTPUT2_sel)-myfx(OUTPUT2_sel[,i])
    }

    
GAP_account[4,2:ncol(GAP_account)]<- Gap_3

setwd('C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS')
write.csv(GAP_account,paste(GAP_name,'.CSV',sep=''),row.names=TRUE)


            #########################################################
            #charging data after second filtering stage if necessary:
            #########################################################

            rawData<- read.table("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/PJC_soil_complete/PJ_2011_soil_complete_07242014.txt", sep = ",", header = TRUE)
            rawData<-rawData[2:17520,]            #only for 2011 data because:all soil_Year_complete datasets begin doy 1 00:30 and finish doy 365 or 366 at 23:30
            rawData <- rawData[, c(1, grep("WC", names(rawData)))]

            OUTPUT<-read.csv("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS/SWC_PJC_2011_Median_filtered_30min.csv",header=T,sep = ",",dec=".")
            OUTPUT<-OUTPUT[,2:ncol(OUTPUT)]

            OUTPUT2<-read.csv("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS/SWC_PJC_2011_REFINED_30min.csv",header=T,sep = ",",dec=".")
            OUTPUT2<-OUTPUT2[,2:ncol(OUTPUT2)]
            OUTPUT2_sel<-OUTPUT2[,6:ncol(OUTPUT2)]

            #GAP_account<-read.csv("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS/SWC_GAPS_PJC_2011_30min.csv",header=T,sep = ",",dec=".")
            GAP_account<-read.csv("C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS/SGAPS_SWC_PJC_2011.csv",header=T,sep = ",",dec=".")
            GAP_account<-GAP_account[,2:ncol(GAP_account)]

            fluxAll<- read.csv("C:/Users/Laura Morillas/Documents/data/PJ_control/Flux_all files/PJ_FLUX_all_2011_Dan.csv",header=T,sep = ",",dec=".")

    

#############################################
#############################################
# PLOTING RESULTS FROM STAGE 2
#   (Refined) OUTPUT2
############################################
############################################

#AUTOMATED LOOP TO SAVE THE diagnostic plots PLOTS
 dev.off()
 
PLOT2_name                        
pdf(file = "C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS/PlotSWC_PJC_2011_REFINED_30min.pdf")
par(mfrow=c(3,1),mar=c(2, 4, 2, 0.5),oma=c(2,0,2,4))

for( i in 2:ncol(rawData)) {
        
plot(timestampDF[,5],rawData[,i],type="l",lwd=2,ylim=c(0,0.5),ylab="SWC",xlab="DOY",main=names(rawData)[i])
lines(timestampDF[,5],OUTPUT[,i],type="l",col="red",lwd=2)
lines(timestampDF[,5],OUTPUT2_sel[,i],type="l",col="green",lwd=2)
par(new=TRUE)                                                                           #This is to create a secondary axis in the same plot
plot(timestampDF[,5],fluxAll[2:nrow(fluxAll),ncol(fluxAll)],type="h",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,max(fluxAll[2:nrow(fluxAll),ncol(fluxAll)],na.rm=T)))        #This is to include precipitation
Yincrease<-(max(fluxAll[2:nrow(fluxAll),ncol(fluxAll)],na.rm=T))/4  
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
Daily_rain<-aggregate(fluxAll[2:nrow(fluxAll),(ncol(fluxAll)-1):ncol(fluxAll)], by=list(fluxAll[2:nrow(fluxAll),ncol(fluxAll)]),FUN=sum, na.rm=TRUE)

OUTPUT1_daily_avg$rain_Tot<-Daily_rain$rain_Tot
OUTPUT1_daily_avg<-OUTPUT1_daily_avg[,2:ncol(OUTPUT1_daily_avg)]

#Saving Daily SWC values
setwd('C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS')
write.csv(OUTPUT1_daily_avg,paste(OUTPUT3_name,'.CSV',sep=''),row.names=TRUE)


#########################
#Accounting the new gaps (daily scale):    
OUTPUT1_daily_avg_sel<-OUTPUT1_daily_avg[,6:ncol(OUTPUT1_daily_avg)]

Gap_4 <- rep(NA, ncol(rawData[,2:ncol(rawData)]))

for( i in grep("WC",names(OUTPUT1_daily_avg_sel))){
     Gap_4[i-1]<- nrow(OUTPUT1_daily_avg_sel)-myfx(OUTPUT1_daily_avg_sel[,i])
    }

    
GAP_account[5,2:ncol(GAP_account)]<- Gap_4

setwd('C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS')
write.csv(GAP_account,paste(GAP_name,'.CSV',sep=''),row.names=TRUE)


########################################
#      Ploting Result 1 THIRD STAGE 
#     (daily averages of refined SWC)
########################################
dev.off()

PLOT3_name

pdf(file = "C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS/PlotSWC_PJC_2011_REFINED_Daily.pdf")
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


INPUT<-OUTPUT1_daily_avg[,c(4,7:ncol(OUTPUT1_daily_avg))]     #A DATASET only including doy and SWC measurements and rain


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
Gap_Rain[i]<-ifelse(Gap_id[i]==1,sum(INPUT[Gap_start2[i]:Gap_end2[i],29],na.rm=T),NA)
                        }
                        
Gap_pos<- (Gap_start[is.na(Gap_start)!=T]+1)
                       
for (i in Gap_pos){

if(Gap_Rain[i]>0)  {SWC_FILLED[Gap_start2[i]:Gap_end2[i]]<-INPUT[Gap_start2[i]:Gap_end2[i],j]}
else {SWC_FILLED[Gap_start2[i]:Gap_end2[i]]<-na.approx(INPUT[Gap_start2[i]:Gap_end2[i],j])}
                        
                       }
                       
SWC_Gapfilled[,j]<-ifelse(Gap_id==1,SWC_FILLED,INPUT[,j])                       
                       
                          }
                          
SWC_Gapfilled$rain_Tot<-INPUT$rain_Tot


setwd('C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS')
write.csv(SWC_Gapfilled,paste(OUTPUT4_name,'.CSV',sep=''),row.names=TRUE)

########################################
#      Ploting Result 2 THIRD STAGE 
#   (daily refined and Gap filled SWC)
########################################

#AUTOMATED LOOP TO SAVE THE diagnostic plots PLOTS
 dev.off()
 
PLOT4_name

pdf(file = "C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS/PlotSWC_PJC_2011_Refined&Gapfilled_Daily.pdf")     #STICK HERE THE NAME OF plot 3
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
    
setwd('C:/Users/Laura Morillas/Documents/data/PJ_control/soil/SWC_FILTERED OUTPUTS')
write.csv(GAP_account,paste(GAP_name,'.CSV',sep=''),row.names=TRUE)


