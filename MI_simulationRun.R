

source("MIEC_simulator.R")
sitelist <- c("AT-Neu","AU-Cpr", "AU-How", "DK-Sor", "FI-Hyy", "FR-Pue", "GF-Guy", "IT-CA1", "US-Los", "US-Ne2")

for (ecd in c(1:10)){
SITE <- sitelist[ecd]
t0 <- Sys.time()

if (SITE=="AT-Neu") {pp <- "FLX_AT-Neu_FLUXNET2015_FULLSET_2002-2012_1-3/FLX_AT-Neu_FLUXNET2015_FULLSET_HH_2002-2012_1-3.csv"; YEAR <- 2010; FREQ <- 1800}
if (SITE=="AU-Cpr") {pp <- "FLX_AU-Cpr_FLUXNET2015_FULLSET_2010-2014_2-3/FLX_AU-Cpr_FLUXNET2015_FULLSET_HH_2010-2014_2-3.csv"; YEAR <- 2012; FREQ <- 1800}
if (SITE=="AU-How") {pp <- "FLX_AU-How_FLUXNET2015_FULLSET_2001-2014_1-3/FLX_AU-How_FLUXNET2015_FULLSET_HH_2001-2014_1-3.csv"; YEAR <- 2011; FREQ <- 1800} 
if (SITE=="DK-Sor") {pp <- "FLX_DK-Sor_FLUXNET2015_FULLSET_1996-2014_2-3/FLX_DK-Sor_FLUXNET2015_FULLSET_HH_1996-2014_2-3.csv"; YEAR <- 2009; FREQ <- 1800}
if (SITE=="FI-Hyy") {pp <- "FLX_FI-Hyy_FLUXNET2015_FULLSET_1996-2014_1-3/FLX_FI-Hyy_FLUXNET2015_FULLSET_HH_1996-2014_1-3.csv"; YEAR <- 2007; FREQ <- 1800}
if (SITE=="FR-Pue") {pp <- "FLX_FR-Pue_FLUXNET2015_FULLSET_2000-2014_2-3/FLX_FR-Pue_FLUXNET2015_FULLSET_HH_2000-2014_2-3.csv"; YEAR <- 2008; FREQ <- 1800}
if (SITE=="GF-Guy") {pp <- "FLX_GF-Guy_FLUXNET2015_FULLSET_2004-2014_2-3/FLX_GF-Guy_FLUXNET2015_FULLSET_HH_2004-2014_2-3.csv"; YEAR <- 2008; FREQ <- 1800}
if (SITE=="IT-CA1") {pp <- "FLX_IT-CA1_FLUXNET2015_FULLSET_2011-2014_2-3/FLX_IT-CA1_FLUXNET2015_FULLSET_HH_2011-2014_2-3.csv"; YEAR <- 2012; FREQ <- 1800}
if (SITE=="US-Los") {pp <- "FLX_US-Los_FLUXNET2015_FULLSET_2000-2014_2-3/FLX_US-Los_FLUXNET2015_FULLSET_HH_2000-2014_2-3.csv"; YEAR <- 2006; FREQ <- 1800}
if (SITE=="US-Ne2") {pp <- "FLX_US-Ne2_FLUXNET2015_FULLSET_2001-2013_1-3/FLX_US-Ne2_FLUXNET2015_FULLSET_HR_2001-2013_1-3.csv"; YEAR <- 2012; FREQ <- 3600}


SCENARIO <- c("S1", "S2", "M5", "L10", "L20")

for (scen in c(1:5)){
	SCEN <- SCENARIO[scen]
	ECimputing_sim(
	pathIN =paste0("~/FLUXNET2015_DATASETS/", pp),
	pathOUT=getwd(),
		site=SITE,
		year=YEAR,
		scenario=SCEN,
		nsim=10,
		M=30,
		SWth = 10)
	}
t1 <- Sys.time()
cat(paste("Elapsed time: ", t1-t0, " Minutes. \n", sep=""))
}

## pathIN: the path where the .csv FLUXNET FULLSET DATA files is stored, e.g.  "C:/FLX_IT-CA1_FLUXNET2015_FULLSET_2011-2014_2-3/FLX_IT-CA1_FLUXNET2015_FULLSET_HH_2011-2014_2-3.csv". Data can be downloaed from https://fluxnet.org
## pathOUT: the path where the results will be stored
## site: character string for the site ID as in FLUXNET database, e.g. IT-CA1
## year: an integer value denoting the calendar year in format YYYY to be analyzed, e.g., 2012
## scenario: one of S1, S2, M5, L10 and L20 as defined in Vitale et al (2018) 
## nsim: an integer value indicating the number of simulations to run for only S1 and S2 scenarios. For M5, L10 and L20 scenarios the number of simulations are 68, 34 and 16, respectively.
## M: an integer value indicating the number of multiple imputations to create for MLR, ADL and PADL MI models.
## SWth: the threshold value (Wm-2) for Potential Radiation used to day-night time regime identification (default 10 Wm-2)

