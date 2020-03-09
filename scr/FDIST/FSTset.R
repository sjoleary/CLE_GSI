FST.set <- function(FST.type){
  if(exists("FSTcalc.sims")){rm(FSTcalc.sims, envir=.GlobalEnv)}
  if(exists("FSTcalc.dataset")){rm(FSTcalc.dataset, envir=.GlobalEnv)}
  if(exists("which.FST")){rm(which.FST, envir=.GlobalEnv)}
  
  if (FST.type=="WC84"){
      FSTcalc.dataset <- function(dataset, diploid,ndig){
        WCtheta.FST.dataset(dataset, diploid, ndig)
      }
  
    FSTcalc.sims <- function(AllCounts){
      WCtheta.FST.Haploids.2Alleles(AllCounts)
    }
    
    which.FST <- function(){"WC84"}    
      
  }else if (FST.type=="CW93"){
      FSTcalc.dataset <- function(dataset, diploid,ndig){
        CW1993.dataset(dataset, diploid, ndig)
      }
  
    FSTcalc.sims <- function(AllCounts){
      CW1993.Beta.Beaumont(AllCounts)
    }  
      which.FST <- function(){"CW93"}  
    
  }else{
    print("Error: Invalid choice for settting FST")
  }   
  myexport(FSTcalc.dataset, FSTcalc.sims, which.FST)
}