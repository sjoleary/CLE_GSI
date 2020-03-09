

LoadCFunctions(sourceWD)

sourceWD <- '~/Dropbox/FODR/'
setwd(sourceWD)
source("Source_FODR_v2.R")
source("Source_FODR_Island.R")

gdata_wd <-"~/Dropbox/FODR/testFODR/"
gdata_name <- "FSTAT_1R_Neut9000_75_20_Set1.dat"
diploid=FALSE
UseMyNullDist=TRUE
totNumDemes <- 75
sample_size <- 20
numberReps <- 1000
MyNullDistFilepath="FSTAT_1R_Neut9000_75_20_Set1.dat_IM_75_20.sim"

#First time simulations done
#a <- IslandOutliers (gdata_name=gdata_name,  diploid=diploid, gdata_wd= gdata_wd, sourceWD=getwd(), totNumDemes= totNumDemes, sample_size=sample_size, UseMyNullDist=FALSE, numberReps=numberReps)

#Use simulations and FST list already done calcHeFST=FALSE, 
#a <- IslandOutliers (gdata_name=gdata_name,  diploid=diploid, gdata_wd= gdata_wd, sourceWD=getwd(), totNumDemes= totNumDemes, sample_size=sample_size, CalcHeFST=FALSE, HeFST_filepath= "FSTAT_1R_Neut9000_75_20_Set1.dat.fst", UseMyNullDist=TRUE, MyNullDistFilepath=MyNullDistFilepath)
#head(a)