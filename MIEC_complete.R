
MIEC_imputing <- function(pathIN, pathOUT, site, year, M, max_resample, SWth){
## pathIN indicates the path where the .csv FLUXNET FULLSET DATA files is stored, e.g.  "C:/FLX_IT-CA1_FLUXNET2015_FULLSET_2011-2014_2-3/FLX_IT-CA1_FLUXNET2015_FULLSET_HH_2011-2014_2-3.csv"
## pathOUT indicates the path where the results will be stored
## site character string for the site ID as in FLUXNET database, e.g. IT-CA1
## year an integer value denoting the calendar year in format YYYY to be analyzed, e.g., 2012
## M an integer value indicating the number of multiple imputations to create for MLR, ADL and PADL MI models.
## max_chain_length
## SWth: the threshold value (Wm-2) for Potential Radiation used to day-night time regime identification (default 10 Wm-2)

library("foreach")
## or backend of choice
library("doMC")
registerDoMC()

 ameliapar <- function(x, m=5, p2s=0, ...) {
	foreach (i = 1:m,
	.combine = "ameliabind",
	.packages="Amelia") %dopar% {
		amelia(x, m=1, p2s=p2s, ...)
		}
	}

library(xts) ## to deal with high-frequency time series
library(forecast) ## to run accuracy function
library(Amelia) ## to run MI
library(strucchange) ## to compute regime detection
library(partitions) ## to build day- and night- time cross-sections
library(data.table) ## to read FLUXNET data
library(REddyProc) ## for MDS algorithm by Reichstein et al (2005)
library(car)
options(scipen = 999) ## to avoid to use scientific notation 
P2S <- 0

cat(paste("Missing Data Imputation for", site, year, "\n"), sep="")
fluxnet_data <- fread(pathIN,  na.strings="-9999", integer64="double", showProgress=FALSE)
timestamp_w <- as.POSIXct(as.character(fluxnet_data[1,1]),format="%Y%m%d%H%M", tz="GMT") + seq(0,FREQ*nrow(fluxnet_data)-1,FREQ)
fluxnet_ts <- xts(fluxnet_data[,-c(1:2)], order.by=timestamp_w)
data2gf <- window(fluxnet_ts, start=as.POSIXct(paste(year,"01010000", sep=""), format="%Y%m%d%H%M", tz="GMT"), end=as.POSIXct(paste(year+1,"01010000", sep=""), format="%Y%m%d%H%M", tz="GMT")-FREQ)
timestamp_s <- time(data2gf)
N <- nrow(data2gf)

dir.create(paste(pathOUT,"/MIFLXres/", site, "_", year, "/complete", sep=""), recursive=TRUE)
output  <- paste(pathOUT,"/MIFLXres/", site, "_", year, "/complete", sep="")
setwd(output)


## VARIABLE SELECTION
Time <- seq(1,N,1)
HH <- round(N/365,0)
HoD <- as.numeric(gl(HH,1,N))
DoY <- as.numeric(gl(N/HH, HH,N))
Rpot <- data2gf$SW_IN_POT
DAYTIME <- which(as.vector(Rpot) > SWth) ## used in accuracy metrics calculation
NIGHTTIME <- which(as.vector(Rpot) <= SWth) ## used in accuracy metrics calculation


if (site != "GF-Guy") ifelse(is.null(data2gf$H_F_MDS), H <- rep(NA,N), H <- as.vector(replace(data2gf$H_F_MDS, which(data2gf$H_F_MDS_QC>0), NA)))

if (site != "GF-Guy") ifelse(is.null(data2gf$LE_F_MDS), LE <- rep(NA,N), LE <- as.vector(replace(data2gf$LE_F_MDS, which(data2gf$LE_F_MDS_QC>0), NA)))

if (site == "GF-Guy") H <- read.table("/Users/mac/Papers/MI-EC_JEI/Review_Analysis/FLX_GF-Guy_FLUXNET2015_FULLSET_2004-2014_2-3/EFDC_L2_Flx_GFGuy_2008_v06_30m.txt", sep=",", header=TRUE, na.strings=-9999)$H

if (site == "GF-Guy") LE <- read.table("/Users/mac/Papers/MI-EC_JEI/Review_Analysis/FLX_GF-Guy_FLUXNET2015_FULLSET_2004-2014_2-3/EFDC_L2_Flx_GFGuy_2008_v06_30m.txt", sep=",", header=TRUE, na.strings=-9999)$LE

ifelse(is.null(data2gf$NEE_VUT_USTAR50), NEE_00 <- rep(NA,N), NEE_00 <- as.vector(replace(data2gf$NEE_VUT_USTAR50, which(data2gf$NEE_VUT_USTAR50_QC>0), NA)))
ind <- which(NEE_00 <= -5)
NEE <- replace(NEE_00, intersect(NIGHTTIME,ind), NA)

ifelse(is.null(data2gf$SW_IN_F_MDS), SWin <- rep(NA,N), SWin <- as.vector(replace(data2gf$SW_IN_F_MDS, which(data2gf$SW_IN_F_MDS_QC>0), NA)))
ifelse(is.null(data2gf$LW_IN_F_MDS), LWin <- rep(NA,N), LWin <- as.vector(replace(data2gf$LW_IN_F_MDS, which(data2gf$LW_IN_F_MDS_QC>0), NA)))
ifelse(is.null(data2gf$NETRAD), Rn <- rep(NA,N), Rn <- as.vector(data2gf$NETRAD))
ifelse(is.null(data2gf$G_F_MDS), G <- rep(NA,N), G <- as.vector(replace(data2gf$G_F_MDS, which(data2gf$G_F_MDS_QC>0), NA)))
ifelse(is.null(data2gf$TA_F_MDS), Ta <- rep(NA,N), Ta <- as.vector(replace(data2gf$TA_F_MDS, which(data2gf$TA_F_MDS_QC>0), NA)))
ifelse(is.null(data2gf$RH), Rh <- rep(NA,N), Rh <- as.vector(data2gf$RH))
ifelse(is.null(data2gf$VPD_F_MDS), VPD <- rep(NA,N), VPD <- as.vector(replace(data2gf$VPD_F_MDS, which(data2gf$VPD_F_MDS_QC>0), NA)))

if (is.null(data2gf$TS_F_MDS_1) & is.null(data2gf$TS_F_MDS_2) & is.null(data2gf$TS_F_MDS_3)) TS <- rep(NA,N)
if (any(!is.null(data2gf$TS_F_MDS_1), !is.null(data2gf$TS_F_MDS_2),!is.null(data2gf$TS_F_MDS_3))) {
	ifelse(is.null(data2gf$TS_F_MDS_1), Ts1 <- rep(NA,N), Ts1 <- as.vector(replace(data2gf$TS_F_MDS_1, which(data2gf$TS_F_MDS_1_QC>0), NA)))
	ifelse(is.null(data2gf$TS_F_MDS_2), Ts2 <- rep(NA,N), Ts2 <- as.vector(replace(data2gf$TS_F_MDS_2, which(data2gf$TS_F_MDS_2_QC>0), NA)))
	ifelse(is.null(data2gf$TS_F_MDS_3), Ts3 <- rep(NA,N), Ts3 <- as.vector(replace(data2gf$TS_F_MDS_3, which(data2gf$TS_F_MDS_3_QC>0), NA)))
	bestTs <- which.min(c(length(which(is.na(Ts1))),length(which(is.na(Ts2))),length(which(is.na(Ts3)))))
	Ts <- as.vector(replace(data2gf[,paste("TS_F_MDS_", bestTs, sep="")], which(data2gf[,paste("TS_F_MDS_",bestTs,"_QC", sep="")]>0), NA))
}

if (is.null(data2gf$SWC_F_MDS_1) & is.null(data2gf$SWC_F_MDS_2) & is.null(data2gf$SWC_F_MDS_3)) SWC <- rep(NA,N)
if (any(!is.null(data2gf$SWC_F_MDS_1), !is.null(data2gf$SWC_F_MDS_2),!is.null(data2gf$SWC_F_MDS_3))) {
	ifelse(is.null(data2gf$SWC_F_MDS_1), SWC1 <- rep(NA,N), SWC1 <- as.vector(replace(data2gf$SWC_F_MDS_1, which(data2gf$SWC_F_MDS_1_QC>0), NA)))
	ifelse(is.null(data2gf$SWC_F_MDS_2), SWC2 <- rep(NA,N), SWC2 <- as.vector(replace(data2gf$SWC_F_MDS_2, which(data2gf$SWC_F_MDS_2_QC>0), NA)))
	ifelse(is.null(data2gf$SWC_F_MDS_3), SWC3 <- rep(NA,N), SWC3 <- as.vector(replace(data2gf$SWC_F_MDS_3, which(data2gf$SWC_F_MDS_3_QC>0), NA)))
	bestSWC <- which.min(c(length(which(is.na(SWC1))),length(which(is.na(SWC2))),length(which(is.na(SWC3)))))
    SWC <- as.vector(replace(data2gf[,paste("SWC_F_MDS_", bestSWC, sep="")], which(data2gf[,paste("SWC_F_MDS_",bestSWC,"_QC", sep="")]>0), NA))
}


ifelse(is.null(data2gf$WS_F), WS <- rep(NA,N), WS <- as.vector(replace(data2gf$WS_F, which(data2gf$WS_F_QC>0), NA)))
ifelse(is.null(data2gf$USTAR), ustar <- rep(NA,N), ustar <- as.vector(data2gf$USTAR))
ifelse(is.null(data2gf$P_F), P <- rep(NA,N), P <- as.vector(data2gf$P_F))
ifelse(is.null(data2gf$TA_ERA), Ta_ERA <- rep(NA,N), Ta_ERA <- as.vector(data2gf$TA_ERA))
ifelse(is.null(data2gf$SW_IN_ERA), SWin_ERA <- rep(NA,N), SWin_ERA <- as.vector(data2gf$SW_IN_ERA))
ifelse(is.null(data2gf$LW_IN_ERA), LWin_ERA <- rep(NA,N), LWin_ERA <- as.vector(data2gf$LW_IN_ERA))
ifelse(is.null(data2gf$VPD_ERA), VPD_ERA <- rep(NA,N), VPD_ERA <- as.vector(data2gf$VPD_ERA))

rm(fluxnet_data)
rm(timestamp_w)
rm(fluxnet_ts)
rm(data2gf)

#############################################################################################################################################################################################################################################################################################################################################################################################################################

## INTERANNUAL REGIMES IDENTIFICATION + FIGURE

#############################################################################################################################################################################################################################################################################################################################################################################################################################

hh <- apply(matrix(NEE,nrow=HH), MARGIN=1, function(x) length(which(is.na(x))))
nee0 <- na.locf(na.locf(na.approx(ts(NEE[seq(which.min(hh),N,HH)], start=1, freq=1), na.rm=FALSE), na.rm=FALSE), fromLast=TRUE, na.rm=FALSE)
nee <- smooth.spline(nee0,nknots=HH/2)$y
bp <- breakpoints(nee~1)
bp0 <-  c(0,(bp$breakpoints-1)*HH,N)
fm1 <- lm(nee~breakfactor(bp, breaks=NULL))
REG <- data.frame(matrix(0, nrow=N, ncol=(length(bp0)-1)))
for (j in 1:(length(bp0)-1)){
	REG[,j] <- replace(REG[,j], seq(bp0[j]+1, bp0[j+1], 1), 1)
	}
colnames(REG) <- paste("Reg", seq(1:(length(bp0)-1)), sep="")
BP_DOY <- format(timestamp_s[bp0[-1]], "%j", tz="GMT")

ifelse(!is.numeric(SWin), SWd <- rep(0,N/HH), SWd <- apply.daily(xts(SWin, order.by=timestamp_s), FUN="mean"))
ifelse(!is.numeric(Ta), Tad <- rep(0,N/HH), Tad <- apply.daily(xts(Ta, order.by=timestamp_s), FUN="mean"))
ifelse(!is.numeric(VPD), VPDd <- rep(0,N/HH), VPDd <- apply.daily(xts(VPD, order.by=timestamp_s), FUN="mean"))
ifelse(!is.numeric(SWC), SWCd <- rep(0,N/HH), SWCd <- apply.daily(xts(SWC, order.by=timestamp_s), FUN="mean"))
ifelse(!is.numeric(P), Pd <- rep(0,N/HH), Pd <- apply.daily(xts(P, order.by=timestamp_s), FUN="sum"))
ifelse(!is.numeric(SWin_ERA), SWERAd <- rep(0,N/HH), SWERAd <- apply.daily(xts(SWin_ERA, order.by=timestamp_s), FUN="mean"))
ifelse(!is.numeric(Ta_ERA), TaERAd <- rep(0,N/HH), TaERAd <- apply.daily(xts(Ta_ERA, order.by=timestamp_s), FUN="mean"))
ifelse(!is.numeric(VPD_ERA), VPDERAd <- rep(0,N/HH), VPDERAd <- apply.daily(xts(VPD_ERA, order.by=timestamp_s), FUN="mean"))

jpeg(paste("Figure_ECORegimes_", site, "_", year, ".jpeg", sep=""),width=480*3, height=480*3)
par(mfrow=c(5,1), mar=c(0,1,0,4), omi=rep(2,4), las=1, cex=1.5, cex.axis=1.5, cex.lab=1.5)
plot(as.vector(NEE), type="l", ylab="", xlab="", xaxt="n", col="grey70")
points(seq(which.min(hh),N,HH), NEE[seq(which.min(hh),N,HH)], pch=21, cex=1.5, bg="grey90")
abline(v=bp0, lty=3, lwd=3)
lines(seq(1,N,HH), fitted(fm1), col = 1, lty=3, lwd=5)

mtext(side=1, at=bp0+HH*10, line=-1.5, c("001", BP_DOY), cex=2)
mtext(side=1, at=bp0[1]+HH*10, line=-2.5, "DoY", cex=2)
mtext(side=2,"NEE", line=5.5, las=0, cex=2)
mtext(side=2,expression((mu*mol*m^-2*s^-1)), line=4, las=0, cex=2)
mtext(side=3, paste(site, year, sep=" "), line=2, cex=3)
mtext(side=3,at=c(bp0[-length(bp0)]+(diff(bp0)/2)), paste("Regime", 1:ncol(REG), sep=" "), cex=2.5)

plot(as.vector(SWd), type="l", xaxt="n", lwd=3, ylab="")
lines(as.vector(SWERAd), col="grey70", lwd=3)
mtext(side=2, "SWin", line=5.5, cex=2, las=0)
mtext(side=2, expression((W*m^-2)), line=4, cex=2, las=0)
abline(v=c(1,as.numeric(BP_DOY)+.5), lty=3, lwd=3)
plot(as.vector(Tad), type="l", xaxt="n", lwd=3, ylab="")
lines(as.vector(TaERAd), col="grey70", lwd=3)
mtext(side=2, "Ta", line=5.5, cex=2, las=0)
mtext(side=2, "(°C)", line=4, cex=2, las=0)
abline(v=c(1,as.numeric(BP_DOY)+.5), lty=3, lwd=3)
plot(as.vector(VPDd), type="l", xaxt="n", lwd=3, ylab="")
lines(as.vector(VPDERAd), col="grey70", lwd=3)
mtext(side=2, "VPD", line=5.5, cex=2, las=0)
mtext(side=2, "(hPa)", line=4, cex=2, las=0)
abline(v=c(1,as.numeric(BP_DOY)+.5), lty=3, lwd=3)
plot(as.vector(SWCd), type="l", xaxt="n", lwd=3, ylim=c(0,max(10,max(as.vector(SWCd),na.rm=TRUE)*1.25)), ylab="")
mtext(side=2, "SWC", line=5.5, cex=2, las=0)
mtext(side=2, "(%)", line=4, cex=2, las=0)
par(new=TRUE)
plot(as.vector(Pd), type="h", col="grey60", lwd=3, xaxt="n", yaxt="n", ylab="")
axis(4)
mtext(side=4, "Prec", line=5, cex=2, las=0)
mtext(side=4, "(mm)", line=3, cex=2, las=0)
abline(v=c(1,as.numeric(BP_DOY)+.5), lty=3, lwd=3)
axis(1, at=c(1,31,60,91,121,152,182,213,244,274,305,335,366), label=c("J","F","M","A","M","J","J","A","S","O","N","D","J"), cex=2)
par(fig = c(0, 1, 0, .9), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE);
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n");
legend("topright", inset=c(0,0),  c("Mean level", paste("Obs at ", format(timestamp_s[which.min(hh)], "%H:%M", tz="GMT"), sep="")), col=1, pch=c(NA,21), pt.bg=c(NA,"grey90"), cex=1.5, lwd=c(7,NA), lty=c(3,NA), bty="n", xpd=TRUE)
legend("topright", inset=c(0,.4),  c("Daily Obs", "Daily ERA"), col=c(1,"grey70"), cex=1.5, lwd=5, bty="n", xpd=TRUE)

dev.off()


#############################################################################################################################################################################################################################################################################################################################################################################################################################

## DAY- NIGHT-TIME REGIMES IDENTIFICATION + FIGURE

#############################################################################################################################################################################################################################################################################################################################################################################################################################

CS <- c()
jpeg(paste("Figure_DiurnalRegimes_", site, "_", year, ".jpeg", sep=""),width=480*3, height=480*4)
par(mfrow=c(ncol(REG),1), mar=c(1,0,2,0), oma=c(6,8,8,4), las=1, cex=1.5, cex.axis=1.5, cex.lab=1.5)
for (j in 1:ncol(REG)){

Rpot_MDC <- apply(matrix(Rpot[which(REG[,j]==1)],nrow=HH), MARGIN=1, function(x) median(x, na.rm=TRUE)) ## simplified version with only Potential Radiation
Daytime <- which(Rpot_MDC > SWth)
Nighttime1 <- seq(1, Daytime[1]-1, 1)
Nighttime2 <- seq(Daytime[length(Daytime)]+1, HH, 1)
D1 <- length(Daytime)
N1 <- length(Nighttime1)
N2 <- length(Nighttime2)
L <- c(N1, D1, N2)

CSname <- c()
for(i in 1:3){
 LL <- L[i]
 if (i==1) {reg="NT1"; asd2 <- LL}
 if (i==3) {reg="NT2"; asd2 <- LL}
 if (i==2) {nCS <- floor(LL/(HH/8)); alpha <- (HH/8); reg="DT"; maxHH <- (HH/6)+1;
 	 ifelse(nCS == 1, asd <- LL, asd <- blockparts(rep(maxHH,nCS),n=LL)); ## f indica il massimo di semiore per CS. asd contiene tutte le combinazioni la cui somma da il numero di semi-ore del daytime e del nighttime
 	 asd1 <- na.omit(t(replace(asd, which(asd<alpha), NA))); ## qui elimino tutte le combinazioni che contengono una valore inferiore ad alpha (il minimo)
 	 ind0 <- apply(asd1, MARGIN=1, function(x) length(which(x==maxHH))); ## qui vedo quali righe contengono il masssimo settato..
 	 ind1 <- which(ind0==0); ## ...e le escludo
 	 asd2 <- asd1[ind1[1],];
 	 if (length(ind1)==0){
 	 	if (dim(table(asd1))==1)  asd2 <- rep(as.numeric(names(table(asd1[1,])))[1],ceiling(as.numeric(table(asd1[1,])[1]/2)))
 	 	if (dim(table(asd1))==2)  asd2 <- c(rep(as.numeric(names(table(asd1[1,])))[1],ceiling(as.numeric(table(asd1[1,])[1]/2))), rep(as.numeric(names(table(asd1[2,])))[2],as.numeric(table(asd1[1,])[2])), rep(as.numeric(names(table(asd1[1,])))[1],as.numeric(table(asd1[1,])[1]/2)))
 	 	}
 	 }

 CSname0 <- paste(reg, c(1:length(asd2)), sep="")
 CSname1 <- rep(CSname0, asd2)
 CSname <- c(CSname, CSname1)
 } 

CS0 <- rep(CSname, length(which(REG[,j]==1))/HH) 
CS <- c(CS, CS0)

xx <- c(rep(which(Rpot_MDC>SWth)[1]-.5,2), rep(which(Rpot_MDC>SWth)[length(which(Rpot_MDC>SWth))]+.5,2))
yy <- c(-100,100,100,-100)
qwe0 <- CS[which(REG[,j]==1)][1:HH]
qwe1 <- rep(NA,47)
for (i in 1:47) {qwe1[i] <- qwe0[i+1]==qwe0[i]}
qwe2 <- which(qwe1==0)+.5
Boxplot(t(matrix(NEE[which(REG[,j]==1)], nrow=HH)), id.n=0, range=3, ylab="", xlab="", xaxt="n", type="n")
abline(h=0, lty=3, lwd=3)
mtext(side=2,"NEE", line=5, las=0, cex=3)
mtext(side=2,expression((mu*mol*m^-2*s^-1)), line=3, las=0, cex=2.5)
mtext(side=3,paste("Regime", j, "from DoY", format(timestamp_s[bp0[j]+1], "%j", tz="GMT"), "to DoY", format(timestamp_s[bp0[j+1]], "%j", tz="GMT"),  sep=" "), line=0.5, las=0, cex=3.5)
axis(1, at=seq(0.5,HH+.5,4), line=0, labels=FALSE, cex=1.5, tick=TRUE, las=2)
axis(3, at=seq(0.5,HH+.5,4), line=0, labels=FALSE, cex=1.5, tick=TRUE, las=2)
polygon(xx,yy, col=rgb(221, 221, 221,300,maxColorValue = 400), border=rgb(221, 221, 221,50,maxColorValue = 255))
#polygon(xx1,yy, col=rgb(221, 221, 221,350,maxColorValue = 400), border=rgb(221, 221, 221,50,maxColorValue = 255))
Boxplot(t(matrix(NEE[which(REG[,j]==1)], nrow=HH)), id.n=0, ylab="", xlab="", xaxt="n", range=3, cex=1.5, lwd=2, add=TRUE)
abline(h=0, lty=3, lwd=3)
abline(v=c(0.5,qwe2,HH+.5),lwd=3)

if(j==1) mtext(side=3, paste(site, year, sep=" "), line=4.5, cex=5)
if(j==ncol(REG)) {
	if(HH==48) axis(1, at=seq(.5,HH+.5,4), line=0, labels=c("00:00","02:00","04:00","06:00","08:00","10:00","12:00","14:00","16:00","18:00","20:00","22:00","24:00"), cex=1.5, tick=FALSE, las=2);
	if(HH==24) axis(1, at=seq(.5,HH+.5,4), line=0, labels=c("00:00","04:00","08:00","12:00","16:00","20:00","24:00"), cex=1.5, tick=FALSE, las=2);
	mtext(side=1, "HoD", cex=3, line=5);
	par(fig = c(.1, .95, 0, .95), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE);
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n");
 	legend("topleft", inset=c(0,0), expression(Rpot>10*W*m^-2), pch=15, col=rgb(221, 221, 221,300,maxColorValue = 400), bty="n", pt.cex=4, cex=2, xpd=TRUE)}
}
dev.off()

#############################################################################################################################################################################################################################################################################################################################################################################################################################

## DATA TRANSFORMING AND PRIORS FOR SOIL PARAMETERS

#############################################################################################################################################################################################################################################################################################################################################################################################################################

PP <- c(); for (i in (HH+1):N) {PP[i] <- sum(P[(i-HH):i])}

ifelse(length(which(is.na(SWC)))/N > .75, SWC_priors <- rnorm(N,0,1), SWC_priors <- na.interp(ts(SWC, start=c(1,1), freq=HH))) ##questa soluzione non è bella ma serve per sopperire ad un baco di Amelia

ifelse(length(which(is.na(Ts)))/N > .75, Ts_priors <- rnorm(N,0,1), Ts_priors <- na.interp(ts(Ts, start=c(1,1), freq=HH)))

#############################################################################################################################################################################################################################################################################################################################################################################################################################

## Set used in MI. "set_stat" refers to those used in the static model specifications while "set_dyn" contain addtional lagged variables

#############################################################################################################################################################################################################################################################################################################################################################################################################################

set_stat <- data.frame(Time, HoD, DoY, CS, H, LE, NEE, SWin, LWin, Rn, G, Ta, Ts, VPD, Rh, SWC, ustar, WS, log1p(PP), Ta_ERA, SWin_ERA, LWin_ERA, VPD_ERA, Ts_priors, SWC_priors)
colnames(set_stat) <- c("Time", "HoD", "DoY", "CS", "H", "LE", "NEE", "SWin", "LWin", "Rn", "G", "Ta", "Ts", "VPD", "Rh", "SWC", "ustar", "WS", "LogP", "Ta_ERA", "SWin_ERA", "LWin_ERA", "VPD_ERA", "Ts_priors", "SWC_priors")
ncol_set_stat <- ncol(set_stat)

set_dyn <- data.frame(Time, HoD, DoY, CS, H, LE, NEE, "H1"=c(NA, H[-N]), "LE1"=c(NA, LE[-N]), "NEE1"=c(NA, NEE[-N]), SWin, LWin, Rn, G, Ta, Ts, VPD, Rh, SWC, ustar, WS,
"SWin1"=c(NA, SWin[-N]),"LWin1"=c(NA, LWin[-N]), "Rn1"=c(NA, Rn[-N]),  "G1"=c(NA, G[-N]), "Ta1"=c(NA, Ta[-N]),"Ts1"=c(NA, Ts[-N]), "VPD1"=c(NA, VPD[-N]), "Rh1"=c(NA, Rh[-N]), "SWC1"=c(NA, SWC[-N]), "ustar1"=c(NA, ustar[-N]), "WS1"=c(NA, WS[-N]), log1p(PP), SWin_ERA, LWin_ERA, Ta_ERA, VPD_ERA, "SWin_ERA1"=c(NA, SWin_ERA[-N]),"LWin_ERA1"=c(NA, LWin_ERA[-N]), "Ta_ERA1"=c(NA, Ta_ERA[-N]), "VPD_ERA1"=c(NA, VPD_ERA[-N]), Ts_priors, SWC_priors)
colnames(set_dyn) <- c("Time", "HoD", "DoY", "CS", "H", "LE", "NEE", "H1", "LE1", "NEE1", "SWin", "LWin", "Rn", "G", "Ta", "Ts", "VPD", "Rh", "SWC", "ustar", "WS", "SWin1","LWin1", "Rn1", "G1", "Ta1","Ts1", "VPD1", "Rh1", "SWC1", "ustar1", "WS1", "LogP", "SWin_ERA", "LWin_ERA", "Ta_ERA", "VPD_ERA", "SWin_ERA1","LWin_ERA1", "Ta_ERA1", "VPD_ERA1", "Ts_priors", "SWC_priors")
ncol_set_dyn <- ncol(set_dyn)


fwrite(data.frame("TIMESTAMP"=format(timestamp_s, "%Y%m%d%H%M", tz="GMT"), set_stat, P, "RPot"=as.vector(Rpot), REG), "DATASET4MI.csv", sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE, append=FALSE)


#############################################################################################################################################################################################################################################################################################################################################################################################################################

## MDS IMPUTATIONS

##############################################################################################################################################################################################################################################################################################################################################################################################################################

cat("Imputing missing values by MDS algorithm: \n")

y <- rep(year, N)
doy <- as.numeric(format(timestamp_s, "%j", tz="GMT"))
if (HH==48) h <- rep(seq(0, 23.5, 0.5), N/HH)
if (HH==24) h <- rep(seq(0, 23, 1), N/HH)
index_qc <- matrix(0, nrow=N, ncol=1)
index0 <- which(is.na(NEE))
index_qc <- replace(index_qc, index0, 2)

set_mds <- data.frame(cbind(y,doy,h, index_qc, NEE, LE, H,  SWin, Ta, Ts, Rh, ustar))
colnames(set_mds) <- c("Year","DoY", "Hour", "qcNEE", "NEE", "LE", "H", "Rg", "Tair", "Tsoil", "rH", "Ustar")
comment(set_mds$Year) <- "--"
comment(set_mds$DoY) <- "--"
comment(set_mds$Hour) <- "--"
comment(set_mds$qcNEE) <- "--"
comment(set_mds$NEE) <- "umolm-2s-1"
comment(set_mds$LE) <- "Wm-2"
comment(set_mds$H) <- "Wm-2"
comment(set_mds$Rg) <- "Wm-2"
comment(set_mds$Tair) <- "degC"
comment(set_mds$Tsoil) <- "degC"
comment(set_mds$rH) <- "%"
comment(set_mds$Ustar) <- "ms-1"

## Export file
write.csv3 <- function(d, file) {
    opts <- options(useFancyQuotes = FALSE)
    on.exit(options(opts))
    h1 <- paste(c(names(d)), collapse = " ")
    h2 <- paste(c(comment(d$Year), comment(d$DoY), comment(d$Hour), comment(d$qcNEE),comment(d$NEE), comment(d$LE),comment(d$H), comment(d$Rg), comment(d$Tair),comment(d$Tsoil), comment(d$rH), comment(d$Ustar)), collapse = " ")
    writeLines(paste(h1, h2, sep = "\n"), file)
    write.table(d, file, sep = " ", append = TRUE, col.names = FALSE, row.names= FALSE, quote=FALSE, eol = "\r\n", na = "-9999")
}
write.csv3(set_mds, "DATASET4MDS.txt")

EddyData.F <- fLoadTXTIntoDataframe('DATASET4MDS.txt')
EddyData.F <- cbind(EddyData.F,VPD=fCalcVPDfromRHandTair(EddyData.F$rH, EddyData.F$Tair))
EddyDataWithPosix.F <- fConvertTimeToPosix(EddyData.F, 'YDH', Year.s="Year", Day.s='DoY', Hour.s='Hour')
EddyProc.C <- sEddyProc$new(site, EddyDataWithPosix.F, c('NEE','LE','H','Rg','Tair','VPD', 'Ustar'), DTS.n=HH)
  
EddyProc.C$sMDSGapFill('Rg', FillAll.b=FALSE) #Fill only the gaps for the meteo condition
EddyProc.C$sMDSGapFill('Tair', FillAll.b=FALSE) #Fill only the gaps for the meteo condition
EddyProc.C$sMDSGapFill('VPD', FillAll.b=FALSE) #Fill only the gaps for the meteo condition
EddyProc.C$sMDSGapFill('NEE', FillAll.b=TRUE) #FillAll.b=TRUE fill all values to estimate flux uncertainties
EddyProc.C$sMDSGapFill('LE', FillAll.b=TRUE) #Fill all values to estimate flux uncertainties  
EddyProc.C$sMDSGapFill('H', FillAll.b=TRUE) #Fill all values to estimate flux uncertainties  
   
FilledEddyData.F <- EddyProc.C$sExportResults()
CombinedData.F <- cbind(EddyData.F, FilledEddyData.F)

fwrite(data.frame(CombinedData.F), "MDSgf.csv", sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)



#############################################################################################################################################################################################################################################################################################################################################################################################################################

## MI IMPUTATIONS

##############################################################################################################################################################################################################################################################################################################################################################################################################################

MI_models <- c("MLR", "ADL", "PADL")

cat("Selected variables and percentage of missing values: \n ")
cat(paste(colnames(set_stat[,5:19]), " ", apply(set_stat[,5:19], MARGIN=2, function(x) round(length(which(is.na(x)))/N*100,1)), "% \n", sep=""))
cat(paste(ncol(REG), " regimes detected. Breakdates at: \n ", sep=""))
cat(paste(format(timestamp_s[c(1,bp0)], "%Y/%m/%d"), "\n", sep=""))


###################################################################################################################################################################################################################################################### MEDIUM LONG GAPS DETECTION
####################################################################################################################################################################################################################################################

LGReg0 <- NULL
NEE1 <- na.locf(na.locf(NEE, maxgap=(HH*5)), maxgap=(HH*5), fromLast=TRUE)

if (length(which(is.na(NEE1)))>0) {
 vec0 <- rep(0,N);
 vec1 <- replace(vec0, which(is.na(NEE1)), 1)
 vec2 <- diff(vec1)
 vec3 <- which(vec2!=0)
 if (vec1[1]==1) vec3 <- c(1,which(vec2!=0))
 if (vec1[N]==1) vec3 <- c(which(vec2!=0),N)
 if (vec1[1]==1 & vec1[N]==1) vec3 <- c(1, which(vec2!=0),N) 
	vec5 <- c()
	vec6 <- c()
	for (i in seq(1,length(vec3),2)){
	vec4 <-  HH*30 ## allungo di 30gg il gap prima e dopo
	vec5 <- c(vec5, ceiling((vec3[i]   - vec4)/HH)*HH+1) ## garantisce che l'inizio del gap sia alla mezzanotte
	vec6 <- c(vec6, ceiling((vec3[i+1] + vec4)/HH)*HH)   ## garantisce che la fine del gap sia alle 23.30
	vec5a <- c(vec5, NA)
	vec6a <- c(NA, vec6)
	ind <- which(vec5a-vec6a<0)
	vec5b <- na.omit(replace(vec5, ind, NA))
	vec6b <- na.omit(replace(vec6, ind-1, NA))
	}
	LGReg0 <- na.omit(cbind(replace(vec5b, which(vec5b<0), 1), replace(vec6b, which(vec6b>N), N)))
	cat(paste(length(vec3)/2, " block(s) of consecutive days missing in NEE time series detected. Shortening regime breakdates and initializing procedure for informative prior estimation. \n", sep=""))
}


## Definizione e adattamento dei REGIMI nel caso di gap > 5 DAYS
if (!is.null(LGReg0)){
  LGREG0 <- matrix(0, nrow=N, ncol=nrow(LGReg0))
  for (j in 1:nrow(LGReg0)){
  LGREG0[,j] <- replace(LGREG0[,j], LGReg0[j,1]:LGReg0[j,2], 1)
  }
  
  LGREG <- matrix(0, nrow=N, ncol=nrow(LGReg0))
  for (j in 1:nrow(LGReg0)){
  LGREG1 <- REG + matrix(rep(LGREG0[,j],ncol(REG)), nrow=N)
  indA <- which.max(apply(LGREG1, MARGIN=2, function(x) length(which(x==2))))
  LGREG[,j] <- replace(LGREG1[,indA], which(LGREG1[,indA]==2),1)
  }

  indB <- which(apply(LGREG, MARGIN=1, FUN="sum")>1)
  if (length(indB)>0){
  for (j in c(1:(ncol(LGREG)-1))){
  indC <- which(apply(LGREG[,c(1:(j+1))], MARGIN=1, FUN="sum")>1)
  LGREG[,j] <- replace(LGREG[,j], indC, 0)
  }
  }
}
 
 RegDummy <- rep(0,N)
 if (length(which(is.na(NEE1)))>0) {
 GAPdummy <- replace(vec0, which(is.na(NEE1)), 1) ## costruisco una dummy con valore 1 per il medium long gap.!!!funziona solo nel caso di 1 medium/long gap. Nel caso di più di un gap non verrano inserite dummy. Questo accadrà per US-Los
 cat(paste(round(length(which(GAPdummy==1))/HH,0), " consecutive days missing in NEE time series detected. \n", sep=""))
 newREG <- REG*GAPdummy
 sumDummy <- apply(newREG, MARGIN=2, FUN="sum")
 vec7 <- which(sumDummy>0)
 vec8 <- intersect(which(REG[,vec7]==1), which(GAPdummy==1))
 if (length(vec8)==length(which(GAPdummy==1))) {RegDummy <- rep(0,N); cat("No Dummy variable for regime shift created. \n")}
 if (length(vec8)<length(which(GAPdummy==1))) {RegDummy <- REG[,vec7[1]]; cat("Dummy variable for regime shift created. \n")}
 }



###################################################################################################################################################################################################################################################### IMPUTATION
####################################################################################################################################################################################################################################################

for (impmod in 1:3){
	ImpMod <- MI_models[impmod]
	file.remove(paste(ImpMod,"overimpNEE.csv", sep="_"))
	file.remove(paste(ImpMod,"overimpLE.csv", sep="_"))
	file.remove(paste(ImpMod,"overimpH.csv", sep="_"))
	file.remove(paste(ImpMod,"Time.csv", sep="_"))
	file.remove(paste(ImpMod,"impNEE.csv", sep="_"))
	file.remove(paste(ImpMod,"impLE.csv", sep="_"))
	file.remove(paste(ImpMod,"impH.csv", sep="_"))
	file.remove(paste(ImpMod,"impALL.csv", sep="_"))
	cat(paste("MI through", ImpMod, "Model \n", sep=" "))
	
	write.table(t(c("Time","NRow","Obs","Fit","Low","Upp","fmv")), paste(ImpMod,"overimpNEE.csv", sep="_"), col.names=FALSE, row.names=FALSE, sep=",", append=FALSE, quote=FALSE)
	write.table(t(c("Time","NRow","Obs","Fit","Low","Upp","fmv")), paste(ImpMod,"overimpLE.csv", sep="_"), col.names=FALSE, row.names=FALSE, sep=",", append=FALSE, quote=FALSE)
	write.table(t(c("Time","NRow","Obs","Fit","Low","Upp","fmv")), paste(ImpMod,"overimpH.csv", sep="_"), col.names=FALSE, row.names=FALSE, sep=",", append=FALSE, quote=FALSE)

###################################################################################################################################################################################################################################################### INFORMATIVE PRIORS IN CASE OF GAP LENGTH > 5 DAYS
####################################################################################################################################################################################################################################################

	if (!is.null(LGReg0)){
		cat(paste("Estimating informative priors for ", ncol(LGREG), " gap(s) > 5 consecutive days in NEE, H and LE time series:\n", sep=""))
		
		set_stat <- data.frame(Time, HoD, DoY, CS, H, LE, NEE, SWin, LWin, Rn, G, Ta, Ts, VPD, Rh, SWC, ustar, WS, log1p(PP), Ta_ERA, SWin_ERA, LWin_ERA, VPD_ERA, Ts_priors, SWC_priors)
		colnames(set_stat) <- c("Time", "HoD", "DoY", "CS", "H", "LE", "NEE", "SWin", "LWin", "Rn", "G", "Ta", "Ts", "VPD", "Rh", "SWC", "ustar", "WS", "LogP", "Ta_ERA", "SWin_ERA", "LWin_ERA", "VPD_ERA", "Ts_priors", "SWC_priors")

		set_dyn <- data.frame(Time, HoD, DoY, CS, H, LE, NEE, "H1"=c(NA, H[-N]), "LE1"=c(NA, LE[-N]), "NEE1"=c(NA, NEE[-N]), SWin, LWin, Rn, G, Ta, Ts, VPD, Rh, SWC, ustar, WS, "SWin1"=c(NA, SWin[-N]),"LWin1"=c(NA, LWin[-N]), "Rn1"=c(NA, Rn[-N]),  "G1"=c(NA, G[-N]), "Ta1"=c(NA, Ta[-N]),"Ts1"=c(NA, Ts[-N]), "VPD1"=c(NA, VPD[-N]), "Rh1"=c(NA, Rh[-N]), "SWC1"=c(NA, SWC[-N]), "ustar1"=c(NA, ustar[-N]), "WS1"=c(NA, WS[-N]), log1p(PP), SWin_ERA, LWin_ERA, Ta_ERA, VPD_ERA, "SWin_ERA1"=c(NA, SWin_ERA[-N]),"LWin_ERA1"=c(NA, LWin_ERA[-N]), "Ta_ERA1"=c(NA, Ta_ERA[-N]), "VPD_ERA1"=c(NA, VPD_ERA[-N]), Ts_priors, SWC_priors)
		colnames(set_dyn) <- c("Time", "HoD", "DoY", "CS", "H", "LE", "NEE", "H1", "LE1", "NEE1", "SWin", "LWin", "Rn", "G", "Ta", "Ts", "VPD", "Rh", "SWC", "ustar", "WS", "SWin1","LWin1", "Rn1", "G1", "Ta1","Ts1", "VPD1", "Rh1", "SWC1", "ustar1", "WS1", "LogP", "SWin_ERA", "LWin_ERA", "Ta_ERA", "VPD_ERA", "SWin_ERA1","LWin_ERA1", "Ta_ERA1", "VPD_ERA1", "Ts_priors", "SWC_priors")

		file.remove("Time.csv")
		file.remove("impNEE.csv")
		file.remove("impLE.csv")
		file.remove("impH.csv")

		for (NLGReg in 1:ncol(LGREG)){
	
	
## DAY- NIGHT-TIME REGIMES IDENTIFICATION FOR MEDIUM LONG GAP TEMPORAL WINDOW
vec0 <- rep(0,N);

Rpot_MDC <- apply(matrix(Rpot[which(LGREG[,NLGReg]==1)],nrow=HH), MARGIN=1, function(x) median(x, na.rm=TRUE)) ## simplified version with only Potential Radiation
Daytime <- which(Rpot_MDC > SWth)#[-1]
Nighttime1 <- seq(1, Daytime[1]-1, 1)
Nighttime2 <- seq(Daytime[length(Daytime)]+1, HH, 1)
D1 <- length(Daytime)
N1 <- length(Nighttime1)
N2 <- length(Nighttime2)
L <- c(N1, D1, N2)

CSname <- c()
for(i in 1:3){
LL <- L[i]
 if (i==1) {reg="NT1"; asd2 <- LL}
 if (i==3) {reg="NT2"; asd2 <- LL}
 if (i==2) {nCS <- floor(LL/6); alpha <- 6; reg="DT"; maxHH <- 9;
 	ifelse(nCS == 1, asd <- LL, asd <- blockparts(rep(maxHH,nCS),n=LL)); ## f indica il massimo di semiore per CS. asd contiene tutte le combinazioni la cui somma da il numero di semi-ore del daytime e del nighttime
	asd1 <- na.omit(t(replace(asd, which(asd<alpha), NA))); ## qui elimino tutte le combinazioni che contengono una valore inferiore ad alpha (il minimo)
	ind0 <- apply(asd1, MARGIN=1, function(x) length(which(x==maxHH))); ## qui vedo quali righe contengono il masssimo settato...
	ind1 <- which(ind0==0); ## ...e le escludo
	asd2 <- asd1[ind1[1],];
	if (length(ind1)==0){
		if (dim(table(asd1))==1)  asd2 <- rep(as.numeric(names(table(asd1[1,])))[1],ceiling(as.numeric(table(asd1[1,])[1]/2)))
		if (dim(table(asd1))==2)  asd2 <- c(rep(as.numeric(names(table(asd1[1,])))[1],ceiling(as.numeric(table(asd1[1,])[1]/2))), rep(as.numeric(names(table(asd1[2,])))[2],as.numeric(table(asd1[1,])[2])), rep(as.numeric(names(table(asd1[1,])))[1],as.numeric(table(asd1[1,])[1]/2)))
		}
 	}
 	
 CSname0 <- paste(reg, c(1:length(asd2)), sep="")
 CSname1 <- rep(CSname0, asd2)
 CSname <- c(CSname, CSname1)
 } 

CS0 <- rep(CSname, length(which(LGREG[,NLGReg]==1))/HH) ## era CS0
##CSlg <- c(CSlg, CS0) qui non è necessario concatenare ma creare un vettore di N elementi con 0 al di fuori del LG e CSlg nel LG
CSlg <- replace(vec0, which(LGREG[,NLGReg]==1), CS0)
	
		
	
## START PRIOR ESTIMATION
if (ImpMod=="MLR") {set0 <- data.frame(set_stat, CSlg, RegDummy); d <- which(Rpot > SWth); n <- which(Rpot <= SWth)}
if (ImpMod=="ADL") {set0 <- data.frame(set_dyn, CSlg, RegDummy); d <- which(Rpot > SWth); n <- which(Rpot <= SWth)}
if (ImpMod=="PADL") {set0 <- data.frame(set_dyn, CSlg, RegDummy); d <- which(substr(CSlg,1,2)=="DT");  n <- which(substr(CSlg,1,2)=="NT")}
ncol_set0 <- ncol(set0)
	
reg_n <- intersect(n, which(LGREG[,NLGReg]==1))
reg_d <- intersect(d, which(LGREG[,NLGReg]==1))
		
			missing_xvar <- round(apply(set0[which(LGREG[,NLGReg]==1),], MARGIN=2, function(x) length(which(is.na(x))))/length(which(LGREG[,NLGReg]==1)),1)
			varzero <- apply(set0[which(LGREG[,NLGReg]==1),], MARGIN=2, function(x) var(x, na.rm=TRUE))
			
			if (ImpMod=="MLR")  {
				idvar_n <- c("Time", "HoD", "DoY", "CS", "SWin", "SWin_ERA", "Ts_priors", "SWC_priors","CSlg")
				idvar_d <- c("Time", "HoD", "DoY", "CS", "Ts_priors", "SWC_priors","CSlg")
				var2exclude <- c(which(varzero[5:ncol_set0]==0)+4, which(missing_xvar[c(5:7)]==1)+5, which(missing_xvar[-c(1:7)]>=.4)+7) ## columns 5-7 refers to EC flux variables and will be excluded only if completely missing. Remaining columns refers to biomet variables and will be excluded only if the percentage of missing is >40%
			}	
				
			if (ImpMod=="ADL") {
				idvar_n <- c("Time", "HoD", "CS", "SWin", "SWin1", "SWin_ERA", "SWin_ERA1", "Ts_priors", "SWC_priors","CSlg")
				idvar_d <- c("Time", "HoD", "CS", "Ts_priors", "SWC_priors","CSlg")
				var2exclude <- c(which(varzero[5:ncol_set0]==0)+4, which(missing_xvar[c(5:10)]==1)+5, which(missing_xvar[-c(1:10)]>=.4)+10) ## columns 5-10 refers to EC flux variables and will be excluded only if completely missing. Remaining columns refers to biomet variables and will be excluded only if the percentage of missing is >40%
			}
							
			if (ImpMod=="PADL") {
				idvar_n <- c("Time", "HoD", "SWin", "SWin1", "SWin_ERA", "SWin_ERA1", "Ts_priors", "SWC_priors", "CS")##qui uso CSlg perchè con CS c'è il rischio di avere CS con numerosità bassa quando prendiamo più regimi
				idvar_d <- c("Time", "HoD", "Ts_priors", "SWC_priors", "CS")
				var2exclude <- c(which(varzero[5:ncol_set0]==0)+4, which(missing_xvar[c(5:10)]==1)+5, which(missing_xvar[-c(1:10)]>=.4)+10) ## columns 5-10 refers to EC flux variables and will be excluded only if completely missing. Remaining columns refers to biomet variables and will be excluded only if the percentage of missing is >40%
			}
			
		
			if (length(var2exclude)==0) set <- data.frame(set0[reg_n,], row.names=NULL)
			if (length(var2exclude)>0)  set <- data.frame(set0[reg_n,-var2exclude], row.names=NULL)

			ts.priors <- rep(NA,4)
			if(any(colnames(set)=="Ts")) {
				missTs <- which(is.na(set[,"Ts"]))
				if(length(missTs)>0) {
					missTs <- which(is.na(set[,"Ts"]));
					ts.m <- set[missTs, "Ts_priors"];
					ts.v <- abs(ts.m)*.05;#sqrt(var(set[,"Ts"], na.rm=TRUE)/10);
					ts.priors <- cbind(missTs, rep(which(colnames(set)=="Ts"), length(missTs)), ts.m, ts.v)
				}
			}	

			SWC.priors <- rep(NA,4)
			if(any(colnames(set)=="SWC")) {
				missSWC <- which(is.na(set[,"SWC"]))
				if(length(missSWC)>0) {
					missSWC <- which(is.na(set[,"SWC"]));
					SWC.m <- set[missSWC, "SWC_priors"];
					SWC.v <- abs(SWC.m)*.05;#sqrt(var(set[,"SWC"], na.rm=TRUE)/10);
					SWC.priors <- cbind(missSWC, rep(which(colnames(set)=="SWC"), length(missSWC)), SWC.m, SWC.v)
				}
			}
			
			my.priors <- na.omit(rbind(ts.priors, SWC.priors))

			sp_n <- 0
			if (length(which(is.na(set[,"NEE"])))/nrow(set) >= .7) ridge_n <- round(0.005*nrow(set),0)
			if (length(which(is.na(set[,"NEE"])))/nrow(set) >= .6 & length(which(is.na(set[,"NEE"])))/nrow(set) < .7) ridge_n <- round(0.0025*nrow(set),0)
			if (length(which(is.na(set[,"NEE"])))/nrow(set) < .6) ridge_n <- NULL

		
			if (dim(my.priors)[1]>0) {
				if (ImpMod=="MLR")  imp_n <- ameliapar(set, idvars=idvar_n, m=M, p2s=P2S, empri=ridge_n, priors=my.priors, incheck=FALSE)
				if (ImpMod=="ADL")  imp_n <- ameliapar(set, idvars=idvar_n, ts="DoY", m=M, p2s=P2S, empri=ridge_n, polytime=sp_n, priors=my.priors, incheck=FALSE)
				if (ImpMod=="PADL") imp_n <- ameliapar(set, idvars=idvar_n, ts="DoY", cs="CSlg", m=M, p2s=P2S, empri=ridge_n, intercs=TRUE, polytime=sp_n, priors=my.priors, incheck=FALSE)
			} 	
			if (dim(my.priors)[1]==0) {
				if (ImpMod=="MLR")  imp_n <- ameliapar(set, idvars=idvar_n, m=M, p2s=P2S, empri=ridge_n, incheck=FALSE)
				if (ImpMod=="ADL")  imp_n <- ameliapar(set, idvars=idvar_n, ts="DoY", m=M, p2s=P2S, empri=ridge_n, polytime=sp_n, incheck=FALSE)
				if (ImpMod=="PADL") imp_n <- ameliapar(set, idvars=idvar_n, ts="DoY", cs="CSlg", m=M, p2s=P2S, empri=ridge_n, intercs=TRUE, polytime=sp_n, incheck=FALSE)
			}	
			cat(paste("Gap ", NLGReg, " during Nighttime: ", imp_n$message, "\n", sep=""))

			if (length(var2exclude)==0) set <- data.frame(set0[reg_d,], row.names=NULL)
			if (length(var2exclude)>0)  set <- data.frame(set0[reg_d,-var2exclude], row.names=NULL)

			ts.priors <- rep(NA,4)
			if(any(colnames(set)=="Ts")) {
				missTs <- which(is.na(set[,"Ts"]))
				if(length(missTs)>0) {
					missTs <- which(is.na(set[,"Ts"]));
					ts.m <- set[missTs, "Ts_priors"];
					ts.v <- abs(ts.m)*.05;#sqrt(var(set[,"Ts"], na.rm=TRUE)/10);
					ts.priors <- cbind(missTs, rep(which(colnames(set)=="Ts"), length(missTs)), ts.m, ts.v)
				}
			}	

			SWC.priors <- rep(NA,4)
			if(any(colnames(set)=="SWC")) {
				missSWC <- which(is.na(set[,"SWC"]))
				if(length(missSWC)>0) {
					missSWC <- which(is.na(set[,"SWC"]));
					SWC.m <- set[missSWC, "SWC_priors"];
					SWC.v <- abs(SWC.m)*.05;#sqrt(var(set[,"SWC"], na.rm=TRUE)/10);
					SWC.priors <- cbind(missSWC, rep(which(colnames(set)=="SWC"), length(missSWC)), SWC.m, SWC.v)
				}
			}
			
			my.priors <- na.omit(rbind(ts.priors, SWC.priors))
			sp_d <- 3
			if (length(which(is.na(set[,"NEE"])))/nrow(set) >= .7) ridge_d <- round(0.005*nrow(set),0)
			if (length(which(is.na(set[,"NEE"])))/nrow(set) >= .6 & length(which(is.na(set[,"NEE"])))/nrow(set) < .7) ridge_d <- round(0.0025*nrow(set),0)
			if (length(which(is.na(set[,"NEE"])))/nrow(set) < .6) ridge_d <- NULL

			if (dim(my.priors)[1]>0) {
				if (ImpMod=="MLR")  imp_d <- ameliapar(set, idvars=idvar_d, m=M, p2s=P2S, empri=ridge_d, priors=my.priors, incheck=FALSE)
				if (ImpMod=="ADL")  imp_d <- ameliapar(set, idvars=idvar_d, ts="DoY", m=M, p2s=P2S, empri=ridge_d, polytime=sp_d, priors=my.priors, incheck=FALSE)
				if (ImpMod=="PADL") imp_d <- ameliapar(set, idvars=idvar_d, ts="DoY", cs="CSlg", m=M, p2s=P2S, empri=ridge_d, intercs=TRUE, polytime=sp_d, priors=my.priors, incheck=FALSE)
			} 	
			if (dim(my.priors)[1]==0) {
				if (ImpMod=="MLR")  imp_d <- ameliapar(set, idvars=idvar_d, m=M, p2s=P2S, empri=ridge_d, incheck=FALSE)
				if (ImpMod=="ADL")  imp_d <- ameliapar(set, idvars=idvar_d, ts="DoY", m=M, p2s=P2S, empri=ridge_d, polytime=sp_d, incheck=FALSE)
				if (ImpMod=="PADL") imp_d <- ameliapar(set, idvars=idvar_d, ts="DoY", cs="CSlg", m=M, p2s=P2S, empri=ridge_d, intercs=TRUE, polytime=sp_d, incheck=FALSE)
			}	
			cat(paste("Gap ", NLGReg, " during Daytime: ", imp_d$message, "\n", sep=""))
			
			set_imp <- array(NA, dim=c(length(which(LGREG[,NLGReg]==1)), ncol(set), M));
			colnames(set_imp) <- colnames(set);
			for (m in (1:M)){
				pre.m <- rbindlist(list(
				imp_n$imputations[[m]],
				imp_d$imputations[[m]]), fill=TRUE, use.names=TRUE)
				set_imp[,,m] <- as.matrix(pre.m[order(pre.m[,"Time"]),])
				rm(pre.m)		
			}

			write.table(set_imp[,"Time",], "Time.csv", append=TRUE, row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)
			write.table(set_imp[,"NEE",], "impNEE.csv", append=TRUE, row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)
			write.table(set_imp[,"LE",], "impLE.csv", append=TRUE, row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)
			write.table(set_imp[,"H",], "impH.csv", append=TRUE, row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)
		}

		time <- read.table("Time.csv", header=FALSE, sep=",")
		impNEE <- read.table("impNEE.csv", header=FALSE, sep=",")
		impLE <- read.table("impLE.csv", header=FALSE, sep=",")
		impH <- read.table("impH.csv", header=FALSE, sep=",")
		TimeLG <- time[order(time[,m]),1]

		NEE_priorsM <- xts(apply(impNEE, MARGIN=1, FUN="mean"), order.by=timestamp_s[TimeLG])
		NEE_priorsSD <- xts(apply(impNEE, MARGIN=1, FUN="sd"), order.by=timestamp_s[TimeLG]) 

		LE_priorsM <- xts(apply(impLE, MARGIN=1, FUN="mean"), order.by=timestamp_s[TimeLG]) 
		LE_priorsSD <- xts(apply(impLE, MARGIN=1, FUN="sd"), order.by=timestamp_s[TimeLG]) 

		H_priorsM <- xts(apply(impH, MARGIN=1, FUN="mean"), order.by=timestamp_s[TimeLG]) 
		H_priorsSD <- xts(apply(impH, MARGIN=1, FUN="sd"), order.by=timestamp_s[TimeLG])

		flux_priors <- merge(Rpot, NEE_priorsM, NEE_priorsSD, LE_priorsM, LE_priorsSD, H_priorsM, H_priorsSD)

		if (ImpMod=="MLR") set_stat <- data.frame(set_stat, coredata(flux_priors[,-1]))
		if (ImpMod=="ADL" | ImpMod=="PADL") set_dyn <- data.frame(set_dyn, coredata(flux_priors[,-1]))
		
		cat("Informative priors estimated successfully. \n \n")
	}
	
	

####### MULTIPLE IMPUTATION

	if (ImpMod=="MLR") {set0 <- set_stat; d <- which(Rpot > SWth); n <- which(Rpot <= SWth)}
	if (ImpMod=="ADL") {set0 <- set_dyn; d <- which(Rpot > SWth); n <- which(Rpot <= SWth)}
	if (ImpMod=="PADL") {set0 <- set_dyn; d <- which(substr(CS,1,2)=="DT");  n <- which(substr(CS,1,2)=="NT")}

	for (NReg in 1:ncol(REG)){
		reg_n <- intersect(n, which(REG[,NReg]==1))
		reg_d <- intersect(d, which(REG[,NReg]==1))
	
		missing_xvar <- round(apply(set0[which(REG[,NReg]==1),], MARGIN=2, function(x) length(which(is.na(x))))/length(which(REG[,NReg]==1)),1)
		varzero <- apply(set0[which(REG[,NReg]==1),], MARGIN=2, function(x) var(x, na.rm=TRUE))
		
		if (ImpMod=="MLR" & is.null(LGReg0))  {
			idvar_n <- c("Time", "HoD", "DoY", "CS", "SWin", "SWin_ERA", "Ts_priors", "SWC_priors")
			idvar_d <- c("Time", "HoD", "DoY", "CS", "Ts_priors", "SWC_priors")
			if (site=="GF-Guy") {idvar_n <- c("Time", "HoD", "DoY", "CS", "SWin", "SWin_ERA", "Ts_priors", "SWC_priors", "Ta_ERA", "VPD_ERA"); idvar_d <- c("Time", "HoD", "DoY", "CS", "Ts_priors", "SWC_priors", "Ta_ERA", "VPD_ERA", "SWin_ERA", "LWin_ERA")}
			
			var2exclude <- c(which(varzero[5:ncol_set_stat]==0)+4, which(missing_xvar==1), which(missing_xvar[8:ncol_set_stat]>=.4 & missing_xvar[8:ncol_set_stat]<1)+7) ## columns 1-7 refers to auxiliary and EC flux variables and will be excluded only if completely missing. Remaining columns refers to biomet variables and will be excluded only if the percentage of missing is >40%
		}
		
		if (ImpMod=="MLR" & !is.null(LGReg0)) {
			ifelse(missing_xvar["NEE_priorsM"]==1, idvar_n <- c("Time", "HoD", "DoY", "CS", "SWin", "SWin_ERA", "Ts_priors", "SWC_priors"), idvar_n <- c("Time", "HoD", "DoY", "CS", "SWin", "SWin_ERA", "Ts_priors", "SWC_priors", "NEE_priorsM", "LE_priorsM", "H_priorsM", "NEE_priorsSD", "LE_priorsSD", "H_priorsSD"))
			ifelse(missing_xvar["NEE_priorsM"]==1, idvar_d <- c("Time", "HoD", "DoY", "CS", "Ts_priors", "SWC_priors"), idvar_d <- c("Time", "HoD", "DoY", "CS", "Ts_priors", "SWC_priors", "NEE_priorsM", "LE_priorsM", "H_priorsM", "NEE_priorsSD", "LE_priorsSD", "H_priorsSD"))
			var2exclude <- c(which(varzero[5:ncol_set_stat]==0)+4, which(missing_xvar==1), which(missing_xvar[8:ncol_set_stat]>=.4 & missing_xvar[8:ncol_set_stat]<1)+7) ## columns 1-7 refers to auxiliary and EC flux variables and will be excluded only if completely missing. Remaining columns refers to biomet variables and will be excluded only if the percentage of missing is >40%
		}
			
		if (ImpMod=="ADL" & is.null(LGReg0)) {
			idvar_n <- c("Time", "HoD", "CS", "SWin", "SWin1", "SWin_ERA", "SWin_ERA1", "Ts_priors", "SWC_priors")
			idvar_d <- c("Time", "HoD", "CS", "Ts_priors", "SWC_priors")
			if (site=="GF-Guy") {idvar_n <- c("Time", "HoD", "CS", "SWin", "SWin1", "SWin_ERA", "SWin_ERA1", "Ts_priors", "SWC_priors", "Ta_ERA", "Ta_ERA1", "VPD_ERA", "VPD_ERA1", "LWin_ERA", "LWin_ERA1"); c("Time", "HoD", "CS", "SWin_ERA","SWin_ERA1", "Ts_priors", "SWC_priors", "Ta_ERA", "Ta_ERA1", "VPD_ERA", "VPD_ERA1", "LWin_ERA", "LWin_ERA1")}
			var2exclude <- c(which(varzero[5:ncol_set_dyn]==0)+4, which(missing_xvar==1), which(missing_xvar[11:ncol_set_dyn]>=.4 & missing_xvar[11:ncol_set_dyn]<1)+10) ## columns 1-10 refers to auxiliary and EC flux variables and will be excluded only if completely missing. Remaining columns refers to biomet variables and will be excluded only if the percentage of missing is >40%
		}
			
		if (ImpMod=="ADL" & !is.null(LGReg0)) {
			ifelse(missing_xvar["NEE_priorsM"]==1, idvar_n <- c("Time", "HoD", "CS", "SWin", "SWin1", "SWin_ERA", "SWin_ERA1", "Ts_priors", "SWC_priors"), idvar_n <- c("Time", "HoD", "CS", "SWin", "SWin1", "SWin_ERA", "SWin_ERA1", "Ts_priors", "SWC_priors", "NEE_priorsM", "LE_priorsM", "H_priorsM", "NEE_priorsSD", "LE_priorsSD", "H_priorsSD"))
			ifelse(missing_xvar["NEE_priorsM"]==1, idvar_d <- c("Time", "HoD", "CS", "Ts_priors", "SWC_priors"), idvar_d <- c("Time", "HoD", "CS", "Ts_priors", "SWC_priors", "NEE_priorsM", "LE_priorsM", "H_priorsM", "NEE_priorsSD", "LE_priorsSD", "H_priorsSD"))
			var2exclude <- c(which(varzero[5:ncol_set_dyn]==0)+4, which(missing_xvar==1), which(missing_xvar[11:ncol_set_dyn]>=.4 & missing_xvar[11:ncol_set_dyn]<1)+10) ## columns 1-10 refers to auxiliary and EC flux variables and will be excluded only if completely missing. Remaining columns refers to biomet variables and will be excluded only if the percentage of missing is >40%
		}

		if (ImpMod=="PADL" & is.null(LGReg0)) {
			idvar_n <- c("Time", "HoD","SWin", "SWin1", "SWin_ERA", "SWin_ERA1", "Ts_priors", "SWC_priors")
			idvar_d <- c("Time", "HoD", "Ts_priors", "SWC_priors")
			if (site=="GF-Guy") {idvar_n <- c("Time", "HoD", "SWin", "SWin1", "SWin_ERA", "SWin_ERA1", "Ts_priors", "SWC_priors", "Ta_ERA", "Ta_ERA1", "VPD_ERA", "VPD_ERA1", "LWin_ERA", "LWin_ERA1"); c("Time", "HoD", "SWin_ERA","SWin_ERA1", "Ts_priors", "SWC_priors", "Ta_ERA", "Ta_ERA1", "VPD_ERA", "VPD_ERA1", "LWin_ERA", "LWin_ERA1")}
			var2exclude <- c(which(varzero[5:ncol_set_dyn]==0)+4, which(missing_xvar==1), which(missing_xvar[11:ncol_set_dyn]>=.4 & missing_xvar[11:ncol_set_dyn]<1)+10) ## columns 1-10 refers to auxiliary and EC flux variables and will be excluded only if completely missing. Remaining columns refers to biomet variables and will be excluded only if the percentage of missing is >40%
		}
			
		if (ImpMod=="PADL" & !is.null(LGReg0)) {
			ifelse(missing_xvar["NEE_priorsM"]==1, idvar_n <- c("Time", "HoD","SWin", "SWin1", "SWin_ERA", "SWin_ERA1", "Ts_priors", "SWC_priors"), idvar_n <- c("Time", "HoD","SWin", "SWin1", "SWin_ERA", "SWin_ERA1", "Ts_priors", "SWC_priors", "NEE_priorsM", "LE_priorsM", "H_priorsM", "NEE_priorsSD", "LE_priorsSD", "H_priorsSD"))
			ifelse(missing_xvar["NEE_priorsM"]==1, idvar_d <- c("Time", "HoD", "Ts_priors", "SWC_priors"), idvar_d <- c("Time", "HoD", "Ts_priors", "SWC_priors", "NEE_priorsM", "LE_priorsM", "H_priorsM", "NEE_priorsSD", "LE_priorsSD", "H_priorsSD"))
			var2exclude <- c(which(varzero[5:ncol_set_dyn]==0)+4, which(missing_xvar==1), which(missing_xvar[11:ncol_set_dyn]>=.4 & missing_xvar[11:ncol_set_dyn]<1)+10) ## columns 1-10 refers to auxiliary and EC flux variables and will be excluded only if completely missing. Remaining columns refers to biomet variables and will be excluded only if the percentage of missing is >40%
		}


		if (length(var2exclude)>0) 	set <- data.frame(set0[reg_n,-var2exclude], row.names=NULL)
		if (length(var2exclude)==0) set <- data.frame(set0[reg_n,], row.names=NULL)

		ts.priors <- rep(NA,4)
			if(any(colnames(set)=="Ts")) {
				missTs <- which(is.na(set[,"Ts"]))
				if(length(missTs)>0) {
					missTs <- which(is.na(set[,"Ts"]));
					ts.m <- set[missTs, "Ts_priors"];
					ts.v <- abs(ts.m)*.05;#sqrt(var(set[,"Ts"], na.rm=TRUE)/10);
					ts.priors <- cbind(missTs, rep(which(colnames(set)=="Ts"), length(missTs)), ts.m, ts.v)
				}
			}	

			SWC.priors <- rep(NA,4)
			if(any(colnames(set)=="SWC")) {
				missSWC <- which(is.na(set[,"SWC"]))
				if(length(missSWC)>0) {
					missSWC <- which(is.na(set[,"SWC"]));
					SWC.m <- set[missSWC, "SWC_priors"];
					SWC.v <- abs(SWC.m)*.05;#sqrt(var(set[,"SWC"], na.rm=TRUE)/10);
					SWC.priors <- cbind(missSWC, rep(which(colnames(set)=="SWC"), length(missSWC)), SWC.m, SWC.v)
				}
			}
			
			my.priors <- na.omit(rbind(ts.priors, SWC.priors))


		if (!is.null(LGReg0)){
			missNEE <- which(is.na(set[,"NEE"]))
			NEE.m <- set[missNEE, "NEE_priorsM"]
			NEE.v <- set[missNEE, "NEE_priorsSD"]
			ifelse(is.null(NEE.m), NEE.priors <- rep(NA,4), NEE.priors <- cbind(missNEE, rep(which(colnames(set)=="NEE"), length(missNEE)), NEE.m, NEE.v))

			missLE <- which(is.na(set[,"LE"]))
			LE.m <- set[missLE, "LE_priorsM"]
			LE.v <- set[missLE, "LE_priorsSD"]
			ifelse(is.null(LE.m), LE.priors <- rep(NA,4), LE.priors <- cbind(missLE, rep(which(colnames(set)=="LE"), length(missLE)), LE.m, LE.v))

			missH <- which(is.na(set[,"H"]))
			H.m <- set[missH, "H_priorsM"]
			H.v <- set[missH, "H_priorsSD"]
			ifelse(is.null(H.m), H.priors <- rep(NA,4), H.priors <- cbind(missH, rep(which(colnames(set)=="H"), length(missH)), H.m, H.v))

			my.priors <- na.omit(rbind(ts.priors, SWC.priors, NEE.priors, LE.priors, H.priors))
		}
 
	        cat(paste(ImpMod,"in Regime", NReg, "during Nighttime.\n", sep=" "))
			sp_n <- 0
			if (length(which(is.na(set[,"NEE"])))/nrow(set) >= .7) {ridge_n <- round(0.005*nrow(set),0);cat(paste("Missing NEE data = ", round(length(which(is.na(set[,"NEE"])))/nrow(set),2)*100, "%. 0.5% of dataset rows added to ridge priors regression.\n", sep=""))}
			if (length(which(is.na(set[,"NEE"])))/nrow(set) < .7) {ridge_n <- round(0.0025*nrow(set),0); cat(paste("Missing NEE data = ", round(length(which(is.na(set[,"NEE"])))/nrow(set),2)*100, "%. 0.25% of dataset rows added to ridge priors regression.\n", sep=""))}
											
			cat("Excluded variable(s): ")
		    if (length(var2exclude)>0) cat(colnames(set0[reg_n,])[var2exclude])
		    if (length(var2exclude)==0) cat("None")
			cat("\nImputing...")

 		ifelse(any(colnames(set)=="NEE"), NEEbound <- t(c(which(colnames(set)=="NEE"), -5, 70)), NEEbound <- t(rep(NA,3)))
		ifelse(any(colnames(set)=="LE"), LEbound <- t(c(which(colnames(set)=="LE"), -100, max(set[,"LE"], na.rm=TRUE)*1.1)), LEbound <- rep(NA,3))
		ifelse(any(colnames(set)=="H"), Hbound <- t(c(which(colnames(set)=="H"), -100, max(set[,"H"], na.rm=TRUE)*1.1)), Hbound <- rep(NA,3))
		ifelse(any(colnames(set)=="SWin"), SWbound <- t(c(which(colnames(set)=="SWin"), 0, max(set[,"SWin"], na.rm=TRUE)*1.1)), SWbound <- rep(NA,3))
		ifelse(any(colnames(set)=="VPD"),VPDbound <- t(c(which(colnames(set)=="VPD"), 0, max(set[,"VPD"], na.rm=TRUE)*1.1)), VPDbound <-  rep(NA,3))
		ifelse(any(colnames(set)=="Rh"), Rhbound <- t(c(which(colnames(set)=="Rh"), 0, 100)), Rhbound <- rep(NA,3))
		ifelse(any(colnames(set)=="WS"), WSbound <- t(c(which(colnames(set)=="WS"), 0, max(set[,"WS"], na.rm=TRUE)*1.1)), WSbound <- rep(NA,3))
		ifelse(any(colnames(set)=="ustar"), USTARbound <- t(c(which(colnames(set)=="ustar"), 0, max(set[,"ustar"], na.rm=TRUE)*1.1)), USTARbound <- rep(NA,3))
		my.bounds <- rbind(NEEbound, LEbound, Hbound, SWbound, VPDbound, Rhbound, WSbound, USTARbound)
	  
			if (dim(my.priors)[1]>0) {
				if (ImpMod=="MLR")	imp_n <- ameliapar(set, idvars=idvar_n, m=M, p2s=P2S, empri=ridge_n, priors=my.priors, incheck=FALSE, bounds=my.bounds, max.resample=max_resample)
				if (ImpMod=="ADL")	imp_n <- ameliapar(set, idvars=idvar_n, ts="DoY", m=M, p2s=P2S, empri=ridge_n, polytime=sp_n, priors=my.priors, incheck=FALSE, bounds=my.bounds, max.resample=max_resample)
				if (ImpMod=="PADL") imp_n <- ameliapar(set, idvars=idvar_n, ts="DoY", m=M, p2s=P2S, empri=ridge_n, cs="CS", intercs=TRUE, polytime=sp_n, priors=my.priors, incheck=FALSE, bounds=my.bounds, max.resample=max_resample)
				} 	

			if (dim(my.priors)[1]==0) {
				if (ImpMod=="MLR")	imp_n <- ameliapar(set, idvars=idvar_n, m=M, p2s=P2S, empri=ridge_n, incheck=FALSE, bounds=my.bounds, max.resample=max_resample)
				if (ImpMod=="ADL")	imp_n <- ameliapar(set, idvars=idvar_n, ts="DoY", m=M, p2s=P2S, empri=ridge_n, polytime=sp_n, incheck=FALSE, bounds=my.bounds, max.resample=max_resample)	
				if (ImpMod=="PADL")	imp_n <- ameliapar(set, idvars=idvar_n, ts="DoY", m=M, p2s=P2S, empri=ridge_n, cs="CS", intercs=TRUE, polytime=sp_n, incheck=FALSE, bounds=my.bounds, max.resample=max_resample)
			}
		cat(paste(imp_n$message, ".\n", sep=""))
		cat("Overimputing...")
	
		oinee_n <- overimpute(imp_n, "NEE")
		oile_n <- overimpute(imp_n, "LE")
		oih_n <- overimpute(imp_n, "H")
		oinee_n_time <- set[,"Time"][oinee_n[,1]]
		oile_n_time <- set[,"Time"][oile_n[,1]]
		oih_n_time <- set[,"Time"][oih_n[,1]]
		fwrite(data.frame(cbind(oinee_n_time, oinee_n)), paste(ImpMod,"overimpNEE.csv",sep="_"), append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
		fwrite(data.frame(cbind(oile_n_time, oile_n)), paste(ImpMod,"overimpLE.csv",sep="_"), append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
		fwrite(data.frame(cbind(oih_n_time, oih_n)), paste(ImpMod,"overimpH.csv",sep="_"), append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
		imp_check0 <- matrix(NA, nrow=M, ncol=3);
		for (i in 1:M){
			imp_check0[i,] <- c(nrow(as.matrix(imp_n$iterHist[[i]])), sum(as.matrix(imp_n$iterHist[[i]][,2])), sum(as.matrix(imp_n$iterHist[[i]][,3])))
		};
		imp_check <- c(min(imp_check0[,1]), max(imp_check0[,1]) , sum(imp_check0[,2]), sum(imp_check0[,3]));
		write.table(t(imp_check), paste(ImpMod,"_impcheck.txt",sep=""), col.names=FALSE, row.names=FALSE, append=TRUE)
		cat("Done. \n \n")


		if(ImpMod=="ADL" | ImpMod =="PADL"){
			NROWn <- nrow(imp_n$imputations$imp1);
			NCOLn <- ncol(set);
			NIGHT_imp <- array(NA, dim=c(NROWn,NCOLn,M));
			colnames(NIGHT_imp) <- colnames(set);
			for (m in (1:M)){
 			NIGHT_imp[,,m] <- as.matrix(imp_n$imputations[[m]][order(imp_n$imputations[[m]][,"Time"]),])
 			};
 			Time.priors0 <- as.numeric(NIGHT_imp[,"Time",1])	;
 			NEE1.priors0 <- apply(NIGHT_imp[,"NEE",], MARGIN=1, function(x) mean(as.numeric(x)));
 			NEE1.priors.xts <- merge(Rpot, xts(NEE1.priors0, order.by=timestamp_s[Time.priors0]))[,-1];
 			rty0 <- replace(set0[,"NEE1"], which(is.na(set0[,"NEE1"])), c(NA,as.vector(NEE1.priors.xts))[which(is.na(set0[,"NEE1"]))]);
 			rty1 <- rty0[reg_d]}
		
		if (length(var2exclude)>0) 	set <- data.frame(set0[reg_d,-var2exclude], row.names=NULL)
		if (length(var2exclude)==0) set <- data.frame(set0[reg_d,], row.names=NULL)
	
		ts.priors <- rep(NA,4)
			if(any(colnames(set)=="Ts")) {
				missTs <- which(is.na(set[,"Ts"]))
				if(length(missTs)>0) {
					missTs <- which(is.na(set[,"Ts"]));
					ts.m <- set[missTs, "Ts_priors"];
					ts.v <- abs(ts.m)*.05;#sqrt(var(set[,"Ts"], na.rm=TRUE)/10);
					ts.priors <- cbind(missTs, rep(which(colnames(set)=="Ts"), length(missTs)), ts.m, ts.v)
				}
			}	

			SWC.priors <- rep(NA,4)
			if(any(colnames(set)=="SWC")) {
				missSWC <- which(is.na(set[,"SWC"]))
				if(length(missSWC)>0) {
					missSWC <- which(is.na(set[,"SWC"]));
					SWC.m <- set[missSWC, "SWC_priors"];
					SWC.v <- abs(SWC.m)*.05;#sqrt(var(set[,"SWC"], na.rm=TRUE)/10);
					SWC.priors <- cbind(missSWC, rep(which(colnames(set)=="SWC"), length(missSWC)), SWC.m, SWC.v)
				}
			}
			
			if (ImpMod=="ADL" | ImpMod =="PADL"){
				NEE1.priors <- rep(NA,4);
				if(any(colnames(set)=="NEE1")) {
				missNEE1 <- which(is.na(set[,"NEE1"]));
				if(length(missNEE1)>0) {
					missNEE1 <- which(is.na(set[,"NEE1"]));
					NEE1.m <- rty1[missNEE1];
					NEE1.v <- abs(NEE1.m)*.05;#sqrt(var(set[,"Ts"], na.rm=TRUE)/10);
					NEE1.priors <- cbind(missNEE1, rep(which(colnames(set)=="NEE1"), length(missNEE1)), NEE1.m, NEE1.v)
				}
			}}

			if(ImpMod=="MLR") {my.priors <- na.omit(rbind(ts.priors, SWC.priors))}
			if(ImpMod=="ADL" | ImpMod =="PADL") {my.priors <- na.omit(rbind(ts.priors, SWC.priors, NEE1.priors))}
			

		if (!is.null(LGReg0)){
			missNEE <- which(is.na(set[,"NEE"]))
			NEE.m <- set[missNEE, "NEE_priorsM"]
			NEE.v <- set[missNEE, "NEE_priorsSD"]
			ifelse(is.null(NEE.m), NEE.priors <- rep(NA,4), NEE.priors <- as.matrix(cbind(missNEE, rep(which(colnames(set)=="NEE"), length(missNEE)), NEE.m, NEE.v)))

			missLE <- which(is.na(set[,"LE"]))
			LE.m <- set[missLE, "LE_priorsM"]
			LE.v <- set[missLE, "LE_priorsSD"]
			ifelse(is.null(LE.m), LE.priors <- rep(NA,4), LE.priors <- cbind(missLE, rep(which(colnames(set)=="LE"), length(missLE)), LE.m, LE.v))

			missH <- which(is.na(set[,"H"]))
			H.m <- set[missH, "H_priorsM"]
			H.v <- set[missH, "H_priorsSD"]
			ifelse(is.null(H.m), H.priors <- rep(NA,4), H.priors <- cbind(missH, rep(which(colnames(set)=="H"), length(missH)), H.m, H.v))

			my.priors <- na.omit(rbind(my.priors, NEE.priors, LE.priors, H.priors))
		}


		    cat(paste(ImpMod, "in Regime", NReg, "during Daytime.\n", sep=" "))
			sp_d <- 3
			if (length(which(is.na(set[,"NEE"])))/nrow(set) >= .7) {ridge_d <- round(0.005*nrow(set),0);cat(paste("Missing NEE data = ", round(length(which(is.na(set[,"NEE"])))/nrow(set),2)*100, "%. 0.5% of dataset rows added to ridge priors regression.\n", sep=""))}
			if (length(which(is.na(set[,"NEE"])))/nrow(set) < .7) {ridge_d <- round(0.0025*nrow(set),0); cat(paste("Missing NEE data = ", round(length(which(is.na(set[,"NEE"])))/nrow(set),2)*100, "%. 0.25% of dataset rows added to ridge priors regression.\n", sep=""))}
							
			cat("Excluded variable(s): ")
		    if (length(var2exclude)>0) cat(colnames(set0[reg_d,])[var2exclude])
		    if (length(var2exclude)==0) cat("None")
			cat("\nImputing...")
			
			#ifelse(any(colnames(set)=="NEE"), NEEbound <- t(c(which(colnames(set)=="NEE"), -100, 70)), NEEbound <- t(rep(NA,3)))
		#ifelse(any(colnames(set)=="LE"), LEbound <- t(c(which(colnames(set)=="LE"), -300, 1000)), LEbound <- rep(NA,3))
		#ifelse(any(colnames(set)=="H"), Hbound <- t(c(which(colnames(set)=="H"), -300, 1000)), Hbound <- rep(NA,3))
		ifelse(any(colnames(set)=="SWin"), SWbound <- t(c(which(colnames(set)=="SWin"), 0, max(set[,"SWin"], na.rm=TRUE)*1.1)), SWbound <- rep(NA,3))
		ifelse(any(colnames(set)=="VPD"),VPDbound <- t(c(which(colnames(set)=="VPD"), 0, max(set[,"VPD"], na.rm=TRUE)*1.1)), VPDbound <-  rep(NA,3))
		ifelse(any(colnames(set)=="Rh"), Rhbound <- t(c(which(colnames(set)=="Rh"), 0, 100)), Rhbound <- rep(NA,3))
		ifelse(any(colnames(set)=="WS"), WSbound <- t(c(which(colnames(set)=="WS"), 0, max(set[,"WS"], na.rm=TRUE)*1.1)), WSbound <- rep(NA,3))
		ifelse(any(colnames(set)=="ustar"), USTARbound <- t(c(which(colnames(set)=="ustar"), 0, max(set[,"ustar"], na.rm=TRUE)*1.1)), USTARbound <- rep(NA,3))
		my.bounds <- rbind(SWbound, VPDbound, Rhbound, WSbound, USTARbound)

			if (dim(my.priors)[1]>0) {
			if (ImpMod=="MLR")  imp_d <- ameliapar(set, idvars=idvar_d, m=M, p2s=P2S, empri=ridge_d, priors=my.priors, incheck=FALSE, bounds=my.bounds, max.resample=max_resample)				
			if (ImpMod=="ADL")  imp_d <- ameliapar(set, idvars=idvar_d, ts="DoY", m=M, p2s=P2S, empri=ridge_d, polytime=sp_d, priors=my.priors, incheck=FALSE, bounds=my.bounds, max.resample=max_resample)
			if (ImpMod=="PADL") imp_d <- ameliapar(set, idvars=idvar_d, ts="DoY", cs="CS", m=M, p2s=P2S, empri=ridge_d, intercs=TRUE, polytime=sp_d, priors=my.priors, incheck=FALSE, bounds=my.bounds, max.resample=max_resample)
			}
			 	
			if (dim(my.priors)[1]==0) {
			if (ImpMod=="MLR")  imp_d <- ameliapar(set, idvars=idvar_d, m=M, p2s=P2S, empri=ridge_d, incheck=FALSE, bounds=my.bounds, max.resample=max_resample)
			if (ImpMod=="ADL")  imp_d <- ameliapar(set, idvars=idvar_d, ts="DoY", m=M, p2s=P2S, empri=ridge_d, polytime=sp_d, incheck=FALSE, bounds=my.bounds, max.resample=max_resample)
			if (ImpMod=="PADL") imp_d <- ameliapar(set, idvars=idvar_d, ts="DoY", cs="CS", m=M, p2s=P2S, empri=ridge_d, intercs=TRUE, polytime=sp_d, incheck=FALSE, bounds=my.bounds, max.resample=max_resample)
		}

	cat(paste(imp_d$message, ".\n", sep=""))
	cat("Overimputing...")

	oinee_d <- overimpute(imp_d, "NEE")
	oile_d <- overimpute(imp_d, "LE")
	oih_d <- overimpute(imp_d, "H")
	oinee_d_time <- set[,"Time"][oinee_d[,1]]
	oile_d_time <- set[,"Time"][oile_d[,1]]
	oih_d_time <- set[,"Time"][oih_d[,1]]
	fwrite(data.frame(cbind(oinee_d_time, oinee_d)), paste(ImpMod,"overimpNEE.csv",sep="_"), append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
	fwrite(data.frame(cbind(oile_d_time, oile_d)),  paste(ImpMod,"overimpLE.csv",sep="_"), append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
	fwrite(data.frame(cbind(oih_d_time, oih_d)),  paste(ImpMod,"overimpH.csv",sep="_"), append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
	imp_check0 <- matrix(NA, nrow=M, ncol=3);
	for (i in 1:M){
		imp_check0[i,] <- c(nrow(as.matrix(imp_d$iterHist[[i]])), sum(as.matrix(imp_d$iterHist[[i]][,2])), sum(as.matrix(imp_d$iterHist[[i]][,3])))
		};
	imp_check <- c(min(imp_check0[,1]), max(imp_check0[,1]) , sum(imp_check0[,2]), sum(imp_check0[,3]));
	write.table(t(imp_check), paste(ImpMod,"_impcheck.txt",sep=""), col.names=FALSE, row.names=FALSE, append=TRUE)
	cat("Done. \n \n")
	
		NROW <- length(which(REG[,NReg]==1))
		NCOL <- ncol(set)
		set_imp <- array(NA, dim=c(NROW,NCOL,M));
		colnames(set_imp) <- colnames(set);
		for (m in (1:M)){
			pre.m <- rbindlist(list(
			imp_n$imputations[[m]],
			imp_d$imputations[[m]]), fill=TRUE, use.names=TRUE)
			set_imp[,,m] <- as.matrix(pre.m[order(pre.m[,"Time"]),])
			rm(pre.m)		
		}
		
		tmp <- as.numeric(set_imp[,"Time",1])
		if (any(colnames(set)=="NEE")) fwrite(data.frame(set_imp[,"NEE",]),  paste(ImpMod,"impNEE.csv",sep="_"), append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
		if (!any(colnames(set)=="NEE")) fwrite(data.frame(matrix(rep(NEE[tmp],M), nrow=NROW, ncol=NCOL)),  paste(ImpMod,"impNEE.csv",sep="_"), append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
		
		if (any(colnames(set)=="LE")) fwrite(data.frame(set_imp[,"LE",]),  paste(ImpMod,"impLE.csv",sep="_"), append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
		if (!any(colnames(set)=="LE")) fwrite(data.frame(matrix(rep(LE[tmp],M), nrow=NROW, ncol=NCOL)),  paste(ImpMod,"impLE.csv",sep="_"), append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
		
		if (any(colnames(set)=="H")) fwrite(data.frame(set_imp[,"H",]),  paste(ImpMod,"impH.csv",sep="_"), append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
		if (!any(colnames(set)=="H")) fwrite(data.frame(matrix(rep(H[tmp],M), nrow=NROW, ncol=NCOL)),  paste(ImpMod,"impH.csv",sep="_"), append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
		

		select_var2exp <- c(which(colnames(set)=="NEE"), which(colnames(set)=="LE"), which(colnames(set)=="H"),  which(colnames(set)=="SWin"), which(colnames(set)=="LWin"), which(colnames(set)=="Rn"), which(colnames(set)=="G"), which(colnames(set)=="Ta"), which(colnames(set)=="Ts"), which(colnames(set)=="VPD"), which(colnames(set)=="Rh"), which(colnames(set)=="SWC"), which(colnames(set)=="ustar"), which(colnames(set)=="WS"))
		if (NReg==1) fwrite(data.frame("TIMESTAMP"=as.numeric(format(timestamp_s[tmp], "%Y%m%d%H%M", tz="GMT")), "RPot"=Rpot[tmp], set_imp[,select_var2exp,1:M],"Prec"=P[tmp]),  paste(ImpMod,"impALL.csv",sep="_"), append=FALSE, row.names=FALSE, col.names=TRUE, quote=FALSE, sep=",")
		if (NReg>1) {tmp_data0 <- fread(paste(ImpMod,"impALL.csv",sep="_"), header=TRUE, sep=",", integer64="numeric"); tmp_data1 <- rbindlist(list(tmp_data0, data.frame("TIMESTAMP"=as.numeric(format(timestamp_s[tmp], "%Y%m%d%H%M", tz="GMT")), "RPot"=Rpot[tmp],set_imp[,select_var2exp,1:M],"Prec"=P[tmp])), fill=TRUE, use.names=TRUE);  fwrite(tmp_data1,  paste(ImpMod,"impALL.csv",sep="_"), append=FALSE, row.names=FALSE, col.names=TRUE, quote=FALSE, sep=",",)}

		} ## END OF MI FOR EACH REGIME
	} ### END OF 3 MI MODELS
} ## END OF IMPUTATION FOR SITE YEAR DATASET

