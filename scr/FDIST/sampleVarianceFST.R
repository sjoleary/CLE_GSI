#code to resample variance due to sample size in FST

### Functions
  DatasetSampleVariance<- function(p.mat, n.mat){
    p.n.mat <-rbind(p.mat,n.mat)
    return(apply(p.n.mat, 2, callMeanSampleVar))
      #the apply function works for columns in the bound p-n matrix
        #each element returned here corresponds to the estimated 
          #sampling variance in FST per locus
      #the columns were bound because the "apply" function only takes one argument
      #the function "callMeanSampleVar" breaks up the bound columns
  }
  
  callMeanSampleVar <- function(p.n.vect){
    n <- length(p.n.vect)/2
    p.obs.vect <- p.n.vect[1:n]
    n.vect <- p.n.vect[(n+1):length(p.n.vect)]
    return(MeanSampleVariance(p.obs.vect, n.vect))
  }  
  
  MeanSampleVariance <- function(p_obs, n_vect){
    # the calculation for a single locus, averaged for 1000 replicates
  	if (length(p_obs)!=length(n_vect)){print ("error in SampleVariance function"); break;}
  	FST_obs <- FSTcalc(p_obs)	
  	return(mean(replicate(1000, SampleVariance(p_obs, n_vect, FST_obs))))
  	}
  	
  SampleVariance <- function(p_obs, n_vect, FST_obs){
  	# The calculation for a single locus and a single replicate	
  	n_pops <- length(p_obs)
  	p_hat <- rbinom(n_pops, n_vect, p_obs)/n_vect
  	FST_hat <- FSTcalc(p_hat)
  	FST_var <- (FST_hat-FST_obs)^2
    return(FST_var)
  	}	
  
    SampleFST <- function(p_obs, n_vect, FST_obs){
    # The calculation for a single locus and a single replicate	
  	n_pops <- length(p_obs)
  	p_hat <- rbinom(n_pops, n_vect, p_obs)/n_vect
  	FST_hat <- FSTcalc(p_hat)
    return(FST_hat)
  	}	
  
    FSTdist <- function(p_obs, n_vect){
      # the calculation for a single locus, averaged for 1000 replicates
      if (length(p_obs)!=length(n_vect)){print ("error in SampleVariance function"); break;}
    	FST_obs <- FSTcalc(p_obs)	
    	return(replicate(1000, SampleFST(p_obs, n_vect, FST_obs)))
  	}
  	
  FSTcalc <- function(p_obs){var(p_obs)/(mean(p_obs)*(1-mean(p_obs)))}	

  
  
# #Example for a single locus
# 	p_obs <- c(0.2, 0.3, 0.5, 0.7, 0.9) #vector of allele frequencies
# 	n_vect <- c(20, 21, 19, 17, 20) #vector of sample sizes for each p calculation
# 	MeanSampleVariance(p_obs, n_vect)
#   
#   #larger sample size
#   n_vect <- c(100,120, 130, 90, 95)
#   MeanSampleVariance(p_obs, n_vect)
# 
# #Example for a dataset with 10 loci in columns and 5 populations in rows
# 	p.obs.mat <- matrix(runif(n=50, 0.1,0.9), ncol=10)
# 	n.mat <- matrix(rpois(50,lambda=20), ncol=10)
# 
#   DatasetSampleVariance(p.obs.mat, n.mat)
#   
#   #confirm that the apply function works for the 1st and 2nd loci
#     MeanSampleVariance(p.obs.mat[,1], n.mat[,1])
#     MeanSampleVariance(p.obs.mat[,2], n.mat[,2])
# 
