

###############################################################################
CW1993.Beta.Eq <- function(AllCounts){
##########################################
### CW93 FST for infinite sample of allele freqs.
### Ref:Cockerham, C. C., and B. S. Weir. 1993. 
###### Estimation of gene flow from F-statistics. Evolution 855-863.
### This function takes a matrix of allele counts, 
### where rows are populations and columns are alleles.
### This function was developed from equations in Cockerham & Weir 1993
### Note that it assumes equal sample sizes
###########################################    
	# AllCounts is a matrix of allele counts at a locus in each population
  # Populations in rows and alleles are in columns
			pops <- rownames(AllCounts)
			numpops <- dim(AllCounts)[1]
			numalleles <- dim(AllCounts)[2]	
			sample.sizes <- rowSums(AllCounts)
      
      if(length(unique(sample.sizes))!=1){
        print("Error in CW1993.Beta.Eq(): sample sizes are not equal")
      }
			
		###### Uncorrected by sample size #####
		#######################################
			a2 <- AllCounts
			p <- AllCounts/sample.sizes
			X <- sum(p^2)
			Y <- sum(colSums(p)^2)
			M <- sample.sizes[1]	
			F0 <- (M*X-numpops)/((M-1)*numpops) 
			F1 <- (Y-X)/(numpops*(numpops-1))
			He <- 1-F1
			FST <- (F0-F1)/(1-F1)
			#FST <- 1 - (1-F0)/(1-F1); an alternative way of writing the previous line
			return(list(numer=F0-F1, He=He, FST=FST))
        # to calculate FST overall loci, sum(numer)/sum(He)
        # because the mean of a ratio is equal to the ratio of the sum(numerator)/sum(denominator)
			}

###############################################################################
CW1993.Beta.Beaumont <- function(AllCounts){
##########################################
### CW93 FST as implemented with sample size correction in FDIST2 
### Ref:Cockerham, C. C., and B. S. Weir. 1993. 
###### Estimation of gene flow from F-statistics. Evolution 855-863.
### Ref: Beaumont, M. A., and R. A. Nichols. 1996. 
###### Evaluating loci for use in the genetic analysis of 
###### population structure. P Roy Soc Lond B Bio 263:1619-1626.
### This function takes a matrix of allele counts, 
### where rows are populations and columns are alleles.
### This function was developed from directly from FDIST2 code
###########################################   
  #AllCounts is a matrix of allele counts at a locus in each population
  #Populations in rows and alleles are in columns

		pops <- rownames(AllCounts)
		numpops <- dim(AllCounts)[1]
		numalleles <- dim(AllCounts)[2]
		sample.sizes <- rowSums(AllCounts)
			
		a2 <- AllCounts
		p <- AllCounts/sample.sizes
		x0 <- sum((rowSums(AllCounts^2)-sample.sizes)/(sample.sizes*(sample.sizes-1)))
				# This is the line that differs from the Beta published in CW1993
				# In CW1993, equal sample sizes are assumed
				# Here, Beaumont applies a sample size correction to each population 
        # sample and then sums over all pops
			
		# The code below is included because it was copied from my.thetacal, 
    # but gives the same x0 as the line above
				#x0<-0
				#for (j in 1:numpops){
				#	x2 <- 0;
				#	for (i in 1:numalleles){
				#		x2<- x2 + AllCounts[j,i]*AllCounts[j,i]
						#print(x2)
				#	}
				#	x0 <- as.numeric(x0 + (x2-sample.sizes[j])/(sample.sizes[j]*(sample.sizes[j]-1)))
				#}
			
			# This code was copied from my.thetacal, 
      # but gives het1 that is the same as (1-F1) in WC1993
			  yy=0
			  for (j in 1:(numpops-1)){
				  for (k in (j+1):numpops){
					  y1=0;
					  for (i in 1:numalleles){
						  y1 <- y1 + AllCounts[j,i]*AllCounts[k,i]
						    #print(y1)
					   }
					  yy <- yy + y1/(sample.sizes[j]*sample.sizes[k])
				  }
			  }
			
			q2<- x0/numpops
			q3 <- 2*yy/(numpops*(numpops-1))
			
			het0 <- 1-q2
			het1 <- 1-q3 #this is what is output as heterozygosity by fdist2
			FST <- (q2-q3)/(1-q3)
  		#FST <- 1-het0/het1 #alternative way of writing    
			return(list(numer=(q2-q3), He=het1, FST=FST))
}

###############################################################################
WCtheta.FST.diploids.book <- function(AllCounts){
##########################################
### WC-Theta FST as implemented (but not used) in FDIST2 
### This function takes a matrix of allele counts, 
### where rows are populations and columns are alleles.
### This function was developed from directly from FDIST2 code
###### Matches "thetacal" function which is included, 
###### but not implemented, with fdist2 code
### nc is the term for the sample size correction
### Estimator taken from page 150 of Weir's book:
### B. S. Weir (1990).  Methods for discrete population genetic data.  
###### Sunderland, MA. Sinauer Publishers.
###########################################   
			pops <- rownames(AllCounts)
			numpops <- dim(AllCounts)[1]
			numalleles <- dim(AllCounts)[2]
			
			sample.sizes <- rowSums(AllCounts)
			a2 <- AllCounts
			p <- AllCounts/sample.sizes
			
		#######################################
			XX <- sum(a2*a2/sample.sizes)
			nbar <- mean(sample.sizes)
			
			q2 <- (XX-numpops)/((nbar-1)*numpops) 
			nc = 1.0/(numpops - 1.0)*(sum(sample.sizes)- sum(sample.sizes^2)/sum(sample.sizes));
			yy <- sum(colSums(a2)^2)
			q3 <- 1/(numpops*(numpops-1)*nbar*nc)*(yy - nbar*(nc-1.0)/(nbar-1.0)*XX) + 
        (nbar-nc)/(nc*(nbar-1.0))*(1.0-1.0/(numpops - 1.0)*XX);
			a0 <- 1-q2
			a1 <- 1-q3
			He <- a1 
			num <- a1-a0
			#FST <- 1 - a0/a1 alternative way of same equation
			FST <- num/He2
			
			return(list(numer=num, He=He, FST=FST))
}

###############################################################################
WCtheta.FST.Haploids.2Alleles<-function(AllCounts){
##########################################
### WC-Theta FST for finite sample of haploid allele freqs
### From page 145-148 of Weir's book
### B. S. Weir (1990).  Methods for discrete population genetic data.  
###### Sunderland, MA. Sinauer Publishers.
### This function is used to calculate FST from the haploid IM simulations
###########################################  
	# Input a matrix of the counts of each allele (columns) in each population (rows)
	if(ncol(AllCounts)!=2){print("Error in WC.FST.Haploids.2Alleles: 
                              only 2 allelic states allowed" )}	
  
	n.pops<-nrow(AllCounts)
	counts1 <- AllCounts[,1]
	sample.sizes <- rowSums(AllCounts)
	n.ave <- mean(as.numeric(sample.sizes))
	n.c = (n.pops*n.ave - sum(sample.sizes^2)/(n.pops*n.ave))/(n.pops-1)	
	p.freqs = counts1/sample.sizes
	p.ave = sum(sample.sizes*p.freqs)/(n.ave*n.pops)
		#note: this differs slightly from mean(p.freqs) in R
	He2 <- 2*p.ave*(1-p.ave)
			
	s2 = sum(sample.sizes*(p.freqs - p.ave)^2)/((n.pops-1)*n.ave)
		#note: this differs slightly from var(p.freqs) in R
  
	T1 <- s2 - 1/(n.ave-1)*(p.ave*(1-p.ave) -(s2*(  n.pops-1)/  n.pops))
	T2 <- (n.c - 1)*(p.ave*(1-p.ave))/(n.ave-1) + 
            (1 + (  n.pops-1)*(n.ave-n.c)/(n.ave-1))*s2/  n.pops

  FST <- T1/T2 
	return(list(numer=T1, denom=T2, He=He2, FST=FST))
}


###############################################################################
WCtheta.FST.dataset <- function(data1, diploid, ndig){
##########################################
### WC-theta FST for finite sample of diploid allele freqs
### return FST and He estimate for each locus in dataset
### Based on Weir and Cockerham's 1984 paper
### Using wc() function from package hierfstat
###########################################

	FST.list <- wc(data1, diploid=diploid)
	FST <- as.numeric(FST.list$per.loc$FST)
  writeLines("Calculating He and FST (Weir & Cockerham 1984) for each locus...")

  numer <- FST.list$sigma.loc$lsiga
  denom <- FST.list$sigma.loc$lsiga + FST.list$sigma.loc$lsigb +
            FST.list$sigma.loc$lsigw
  
	FST.overall <- as.numeric(FST.list$FST)
  FIS.overall <- as.numeric(FST.list$FIS)
  He <- as.numeric(getHe.dataset(data1, diploid, ndig))
	locnames <- colnames(data1)[-1]

	return(list(FST.overall=FST.overall, FIS.overall=FIS.overall, 
              FST.emp =data.frame(locnames =locnames, He=He, FST=FST, 
              numer=numer, denom=denom)))
}

###############################################################################
mean.of.ratios <- function(numer, denom){
  #numer and denom are each vectors
  sum(numer)/sum(denom)
}

###############################################################################
CW1993.dataset <- function(data1, diploid, ndig){
########################   
### getCW1993 from a dataframe in FSTAT format
### df is the dataframe
### diploid is TRUE for diploids and FALSE for haploids
### ndig is the number of digits used to code genotypes
########################
  pops <- data1[,1]
  nloc <- ncol(data1)-1
	locnames <- colnames(data1)[-1]

	F0F1 <- matrix( unlist(apply(data1[,-1], MARGIN=2,get.CW93.locus, diploid,
                               ndig, pops)), 
                  ncol= 3, byrow=TRUE)
	colnames(F0F1) <- c("Numer", "He", "Fst")
  F0F1[which(is.na(F0F1))]<-NA
  numer<-F0F1 [,1] 
  He <- F0F1[,2]
  FST <- F0F1[,3]
  FST.overall <- mean.of.ratios(numer,He)

  return(list(FST.overall=FST.overall, FST.emp =data.frame( He=He, FST=FST, 
              numer=numer, denom=denom, locnames =locnames)))

}

get.CW93.locus <- function(vector, diploid, ndig, pops){
### Internal function to CW1993.dataset
    	countMat <- getAllCounts(pops, vector, diploid, ndig)
      return(CW1993.Beta.Beaumont(countMat))
}      