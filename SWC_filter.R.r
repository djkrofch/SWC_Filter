################################################################################
# swcFilter.v02.R edited Jul 22 2014
################################################################################
#
# USAGE: This is not a production script nor should it be considered free from
# 	 faults or errors. For testing purposes only.
#
#	 This script quickly assesses the capability of a wavelet filter to
#	 reduce the noise and smooth the time series of the SWC data collected
#	 by the CS616 probes at the PJ sites. These data have been characterized
#	 as having a consistent 'electrical noise' signal which results in a 
#	 seemingly random downward spike in the data at intermittent intervals
#	 but apparently consistent magnitudes. Given the desire to preserve the
# 	 perceived 'real' positive spikes in the signal due to precipitation,
#	 but not negatively affect the mean signal during dry, noisey periods,
#	 a wavelet filter may produce better results than a median, variance,
#	 or mean moving average window (time domain) based approach. 
#
#
# For inquiries about this script's use or function, contact the authors at:
#	Dan Krofcheck      krofcheck@gmail.com
#
#################################################################################
#################################################################################

### Depends package Reshape

#install.packages("ggplot2")
#install.packages("reshape")
#install.packages("robfilter")

library(ggplot2)
library(reshape)
library(zoo)
library(robfilter)

myfx<- function (x){ sum(!is.na(x))} #function to account the number of values that are not NA

### Read in an example CSV containing raw SWC data from PJC
                        
testData <- read.table("C:/Users/Laura Morillas/Documents/data/PJ_control/PJC_soil_complete/PJ_2011_soil_complete_07242014.txt", sep = ",", header = TRUE)
                       
testData<-testData[2:17520,]            #only for 2011 data because:all soil_Year_complete datasets begin doy 1 00:30 and finish doy 365 or 366 at 23:30
fluxAll  <- read.csv("C:/Users/Laura Morillas/Documents/data/PJ_control/Flux_all files/PJ_FLUX_all_2011_Dan.csv",header=T,sep = ",",dec=".")
                                        #all Flux_all files begin doy 1 00:00 and finish doy 365 or 366 at 23:30
                                        
### Select only the columns that contain SWC content. Searching here for *WC* as output
### by the CR23x. Also include the matlab datenum timestamp for POSIX conversion
testData <- testData[, c(1, grep("WC", names(testData)))] #grep functions look for the colnames that includes "WC"

### Include precipitation in the working dataset (from flux all files): 
#testData$PRECIP<-fluxAll[2:nrow(fluxAll),ncol(fluxAll)]   #This will work always because soil_Year_complete datasets begin doy 1 00:30 and finish doy 365 or 366 at 23:30 
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

timestamp <- matlab2POS(testData$timestamps)

timestampDF <- data.frame(ts = timestamp, year =  timestamp$year + 1900, month = timestamp$mon,
	mday = timestamp$mday, doy = (timestamp$yday)+1, 
	hrmin = (timestamp$hour + timestamp$min / 60))
	
#Creating a matrix to record output files 

OUTPUT0<- mat.or.vec(nrow(timestampDF), nc=length(names(testData)))
OUTPUT<-data.frame(OUTPUT0)
names(OUTPUT)<-names(testData[,c(1:dim(testData)[2])])
OUTPUT[,1]<-testData$timestamp

GAP_account0<- mat.or.vec(nr=3, nc=length(names(testData)))
GAP_account<-data.frame(GAP_account0)
names(GAP_account)<- c("Gap_level",names(testData[,c(2:dim(testData)[2])]))
GAP_account[,1]<-c("Gap_0","Gap_1","Gap_2")



########################
### CURRENT LOOP START
for(j in 1:ncol(testData)){

colname <- names(testData)[j+1]           #To determine the probe that we are going to check
currentCol <- data.frame(testData[, j+1])
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
                                                                #AQUI METER UNA COLUMNA QC PARA CADA SENSOR CON UN DIFERENTE VALOR SEGUN EL TIPO DE GAP FILLING. eSTE PRIMER PASO SERIA QC=1
                                                                
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

write.csv(OUTPUT,paste('PJC_11_FilteredSWC_WS50','.CSV',sep=','),row.names=TRUE)
write.csv(GAP_account,paste('PJC_11_SWC_GapFilteraccounter_WS50','.CSV',sep=','),row.names=TRUE)


#######################
#Ploting results:
#######################


#AUTOMATED LOOP TO SAVE THE diagnostic plots PLOTS
 dev.off()
                        
#jpeg(file = "C:/Users/Laura Morillas/Documents/R/SWC filtering/Diagnostic plots/Testplotfinal.jpg",width=1500,height=1000)
#png(file = "C:/Users/Laura Morillas/Documents/R/SWC filtering/Diagnostic plots/Testplotfinal.png") 
#THIS TWO, JPEG AND PNG, TURN CRAZY WITH THE LOOP AND RECORD ONLY SOME PROBES

pdf(file = "C:/Users/Laura Morillas/Documents/R/SWC filtering/Diagnostic plots/prueba.pdf")
par(mfrow=c(3,1),mar=c(2, 4, 2, 0.5),oma=c(2,0,2,4))

for( i in 2:ncol(testData)) {
        
plot(timestampDF[,5],testData[,i],type="l",lwd=2,ylim=c(0,0.5),ylab="SWC",xlab="DOY",main=names(testData)[i])
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
#Paying with medianm filter window size:
############################################

OUTPUT_WS30<-OUTPUT
OUTPUT_WS50<-OUTPUT
OUTPUT_WS70<-OUTPUT
OUTPUT_WS100<-OUTPUT
OUTPUT_WS150<-OUTPUT


#To choose specific time periods for the plot:
RA<-c(1:nrow(timestampDF))
timestampDF_col<-cbind(timestampDF,RA)
TE<-range(timestampDF_col[timestampDF_col$doy>=200 & timestampDF_col$doy<=270,7])

Tsubset<-c(TE[1]:TE[2])

pdf(file = "C:/Users/Laura Morillas/Documents/R/SWC filtering/Diagnostic plots/PJC_2011_SWCfilter_WScompDetail30to70_2.pdf")
par(mfrow=c(3,1),mar=c(2, 4, 2, 0.5),oma=c(2,0,2,4))

for( i in 2:ncol(testData)) {
        
plot(timestampDF_col[Tsubset,7],testData[Tsubset,i],pch=8,type="p",lwd=2,ylim=c(0,0.5),ylab="SWC",main=names(testData)[i]) #xlab="DOY",
lines(timestampDF_col[Tsubset,7],OUTPUT_WS30[Tsubset,i],pch=4,cex=0.3 ,type="p",col="pink",lwd=2)
lines(timestampDF_col[Tsubset,7],OUTPUT_WS50[Tsubset,i],pch=4,cex=0.3 ,type="p",col="green",lwd=2)
lines(timestampDF_col[Tsubset,7],OUTPUT_WS70[Tsubset,i],pch=4,cex=0.3 ,type="p",col="purple",lwd=2)
lines(timestampDF_col[Tsubset,7],OUTPUT_WS100[Tsubset,i],pch=4,cex=0.3 ,type="p",col="orange",lwd=2)
lines(timestampDF_col[Tsubset,7],OUTPUT_WS150[Tsubset,i],pch=4,cex=0.3 ,type="p",col="red",lwd=2)
par(new=TRUE)                                                                           #This is to create a secondary axis in the same plot
plot(timestampDF_col[Tsubset,7],fluxAll[Tsubset,ncol(fluxAll)],type="h",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,max(fluxAll[2:nrow(fluxAll),ncol(fluxAll)],na.rm=T)))        #This is to include precipitation
Yincrease<-(max(fluxAll[2:nrow(fluxAll),ncol(fluxAll)],na.rm=T))/4  
axis( side=4)
mtext("Precip",side=4,line=3,cex=0.8)                                                                           #4 means that the rain axis is going to be writen in the right hand of the plot
legend("topright",col=c("black","pink","green","purple","orange","red","blue"),lty=1,legend=c("raw SWC","Filtered SWC WS30","Filtered SWC WS50"
,"Filtered SWC WS70","Filtered SWC WS100","Filtered SWC WS150","Precipitation"))

     } 
     
dev.off()
dev.off()
dev.off()

#Note a value >150 for the window size parameter in the mediam filter reduces too much the higher values of SWC as a result of rain.
# window size of 30 or 50 seems to be the better option allowing the SWC keep realistic under fast SWC increases as a result to rain




######################################################
#####################################################
### To built a noisy quality acounter after filter SWC (ws=50)
#####################################################
#####################################################


#TESTING WITH ONLY ONE COLUMN:
df<-OUTPUT[,c(1,6)]
df$PRECIP<-fluxAll[2:nrow(fluxAll),ncol(fluxAll)]
 
windowsize = 47                           #1 DAY WINDOW
var_df<-rep(NA,nrow(df))

Q<-round(windowsize)                  #THE FIRST ROW INCLUDED IN WINDOW
Z<-(round(nrow(df)-windowsize/2)-1)   #THE LAST ROW INCLUDED IN WINDOW

#tRYING WITH STANDARD DEVIATION

for (i in Q:Z){

window <- df[(i - windowsize):i, ]
windowP <- df[(i - windowsize/2):(i+2),]

if (sum(windowP[,3],na.rm=T) > 0 ){   #& !is.na(sum(window[,3]) { 
      var_df[i]<-0 
      }
else {      
var_df[i] <- sd(window[,2]) #use var() or standard error and pay with the best
}

}

var_df_2<-ifelse( var_df >=(0.2*max(var_df,na.rm=T)),1,0)


#tO CHECK quality acounter results

par(mfrow=c(2,1))

Tsubset<-c(9600:13007)  #For monsoon
Tsubset<-c(10000:11000) #For monsoon MORE DETAIL
Tsubset<-c(1:nrow(df))  #For all data

plot(df[Tsubset,1],var_df_2[Tsubset],col="green",type="h",ylim=c(0,max(var_df,na.rm=T)))#,ylim=c(0,0.005)
lines(df[Tsubset,1],var_df[Tsubset],col="orange",type="h")#,ylim=c(0,0.005)
par(new=TRUE)
plot(df[Tsubset,1],df[Tsubset,2],main=names(df)[2])
axis(4)
par(new=TRUE)
plot(df[Tsubset,1],fluxAll[Tsubset,ncol(fluxAll)],type="h",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,max(fluxAll[2:nrow(fluxAll),ncol(fluxAll)],na.rm=T)))        #This is to include precipitation


#par(new=TRUE)
#plot(df[Tsubset,1],var_df[Tsubset],col="green",type="h")#,ylim=c(0,0.005)



#####
#To generate a quality column for each probe   (under construction yet)
####
OUTPUT_QC_0<- mat.or.vec(nrow(OUTPUT), nc=length(names(OUTPUT)))
OUTPUT_QC<-data.frame(OUTPUT_QC_0)
names(OUTPUT)<-names(testData[,c(1:dim(testData)[2])])
OUTPUT[,1]<-testData$timestamp






###Dan
### Now we have filled all the small gaps using linear approximation. The more complicated gaps
### require more robust fitting. Lets first find all the gaps and plot them.
gaps   <- data.frame()
buffer <- 5 
count  <- 1
for(i in 1:nrow(timestampDF)){
	if(is.na(SWCdf_2$SWC_fill2[i])){
		gapstart <- (i-1)
		while(!is.na(SWCdf_2$SWC_fill2[i])){
			print(i)
			print(count)
			print(gapstart)
			i = i + 1
		}
		gapstop <- (i) 
		gaps <- rbind(gaps, cbind(SWCdf_2$SWC_fill2[(gapstart - buffer) : (gapstop + buffer)], 
			SWCdf_2$ts[(gapstart - buffer) : (gapstop + buffer)], count))

		count <- count + 1
		print(i)
	}
}

names(gaps) <- c('SWC','ts','count')


layout(matrix(c(1,1,1,1,1,1),2,3, byrow = TRUE))
for(i in max(gaps$count)){
	
}
