
get.pval <- function(FST.obs, He.obs, FST.sim, He.sim, He.bin=0.04){
########################   
##### Get p value for a single observed FST, based on simulated values
########################
	# FST.obs is a single value for the observed locus
	# He.obs is a single value for the observed locus
	# FST.sim is the vector of simulated FSTs
	# He.sim is the vector of simulated He values
	# He.bin is the bin around the observed He for which the p-value will be estimated
    # default is 0.04 following Excoffier et al (2009)
	
### Get the list of FST values that are within the Heterozygosity-bin of the observed He
	
	MinHe <- min(He.sim)
	MaxHe <- max(He.sim)
	if ((He.obs-MinHe)<He.bin/2){ 
    #if the observed He is near the lower bound, then arrange the bin
			He.ind <- which(He.sim<(MinHe+He.bin)  & He.sim>MinHe)
	}else if((MaxHe-He.obs)<He.bin/2){
    #if the observed He is near the upper bound, then arrange the bin
			He.ind <- which(He.sim<(MaxHe)  & He.sim>(MaxHe-He.bin))
	}else{
			He.ind <- which(He.sim<(He.obs+He.bin/2) & He.sim>(He.obs-He.bin/2))
	}	
	
	FST.He <- FST.sim[He.ind]	
		
### Fit the simulated FSTs to a Johnson distribution and get the density
	minFST <- min(c(min(FST.sim), FST.obs))
	FST <- seq(minFST,1.01,by=0.0001)
		# just to make sure "1" is included in sequence
	FST.dens <- dJohnson(FST, parms=JohnsonFit(FST.He))
	
	FST.dens[which(FST.dens=="NaN")]=0
		# 0 probabilities are undefined in the function

### Normalize the density and get the p-value for the observed FST	
	FST.dens2 <- FST.dens/sum(FST.dens)
	FST.ind <- max(which(FST<=FST.obs))
	p.val<- sum(FST.dens2[1:FST.ind])
	return(p.val)
}


get.pval.dataset <- function(FST.empDF, FST.simDF, He.bin=0.04, 
                             write.progress=TRUE){
##########################################
### Applies the get.pval function to all loci in the dataset
### FST.empDF is a dataframe with a column for the observed He and FST for each locus
    # it is OK for this dataframe to have other columns, P-values will just be appended
### FST.simDF is a dataframe with a column for simulated He and FST for several thousand loci
###### note that these two dataframes should not have the same number of rows
### write.progress will write a line every 1000 loci
###########################################  
 
	FST.empVect <- FST.empDF$FST
	He.empVect <- FST.empDF$He
	
	FST.simVect <- FST.simDF$FST
	He.simVect <- FST.simDF$He
	
	ntimes <- length(FST.empVect)
	p.val <- rep(NA, ntimes)
	for (i in 1:ntimes){
			p.val[i] <- get.pval(FST.obs=FST.empVect[i], He.obs=He.empVect[i],  
                           FST.sim=FST.simVect, He.sim=He.simVect, He.bin=He.bin)
    if (write.progress=="TRUE"){  
  		if (i%%1000==0) {writeLines(paste("P-values calc'd for ", i, "of", 
                                        ntimes, "loci", sep=" "))
  		}
    }
	}

	out1 <- as.data.frame(cbind(FST.empDF, p.val.cum=p.val))
	return(out1)
}


correct.pval.dataframe <- function(dataframe, p.colName="p.val.cum", p.colNum){
##########################################
### Correct a list of cumulative p-values to indicate tails,
###### and indicate significance at the Bonferroni, FDR=0.05, and FDR=0.01 levels
###### Depends on Storey's qvalue.R function
### infilepath is path to the file, if it is not yet loaded as a dataframe
### assumes the column in the dataframe is named p.val.cum
    ### if the column has another name, please specify name and number of the column 
    ### the p-value column should be based on cumulative probabilities 
    ###	(i.e. starting at 0 in the left tail and ending at 1 in the right tail)
### In the output, "L" indicates left-tail, and "R" indicates right-tail	
### if write.outfile==TRUE
    ### will write to a file in outfilepath, but with the extension ".Cpval"
    ### if outfilepath is NOT specified, will write to infilepath with the extension ".Cpval"
    ### if neither outfile path or infilepath is specified, will give an error

  if (p.colName!="p.val.cum"){
    p.val <- assign(p.colName, dataframe[,p.colNum])
  }else{
    p.val <- dataframe$p.val.cum
  }
		
	L.p <- p.val
	R.p <- 1-L.p
	num.obs <- length(L.p)
	Bonf.p <- 0.05/num.obs
	
  ### Get qvalues for left hand side of distribution
  	q2 <- qvalue(L.p)
  	L.q <- q2$qvalues
  		
  ### Get qvalues for right hand side of distribution	
  	q3 <- qvalue(R.p)
  	R.q <- q3$qvalues
  			
  ### Convert p-values and q-values to the right and left sides of the distribution
  	p.val.tail <- L.p
  	p.val.tail[L.p>0.5] <- R.p[L.p>0.5] 
  	
  	qval <- L.q
  	maxq <- min(q2$pi0, q3$pi0) 
      #if pi0 is different for the left and right sides, take the minimum
  	qval[L.q<maxq] <- L.q[L.q<maxq]
  	qval[R.q<maxq] <- R.q[R.q<maxq]
  		
  	Tail <- L.p
  	Tail[L.p>0.5] <- "R"
  	Tail[L.p<=0.5] <- "L"
  	tail <- as.factor(Tail)
  		
  	Bonf <- rep(FALSE, num.obs)
  	Bonf[p.val.tail<Bonf.p] <- TRUE
  
  	FDR.01 <- rep(FALSE, num.obs)
  	FDR.01[qval<0.01] <-  TRUE
  	
  	FDR.05 <- rep(FALSE, num.obs) 
  	FDR.05[qval<0.05] <-  TRUE
  	
  	out <- data.frame(dataframe, tail=tail, p.val.tail=p.val.tail, qval.tail=qval, 
                      Bonf=Bonf, FDR.01=FDR.01, FDR.05=FDR.05)

	return(out)
}

write.all.PvalueInfo.dataset <- function(FST.emp, FST.sim, outfilepath, 
                                       write.outfile=TRUE){
    if (write.outfile==TRUE & length(outfilepath)==0){
      writeLines("Error in write.all.PvalueInfo.dataset.  No outfilename specified."); break
    }
  if (write.outfile==TRUE & file.exists(outfilepath)){
    writeLines("Error in write.all.PvalueInfo.dataset.  Outfilename already exists."); break
  } 
    
   p.vals<- get.pval.dataset(FST.emp, FST.sim)
   p.vals.C <- correct.pval.dataframe(p.vals)  

   write.table(p.vals.C,outfilepath, row.names=FALSE)
   return(p.vals.C)
}

ReturnCorrectedPVal.BS <- function(filepath){
### BayeScan outputs a file with q-values, but not p-values, already corrected for  
### this function returns corrected p values for a Bayescan outfile

  bs <- read.table(filepath, header=TRUE)
		p.valC <- as.numeric(1-bs$prob)
		num.obs <- length(bs$prob)
		cent <- mean(as.numeric(bs$fst))
	
		Tail <- rep(NA, num.obs)
		Tail[bs$fst<cent] <- "L"
		Tail[bs$fst>cent] <- "R"
		Tail <- as.factor(Tail)
		
		Bonf.CO <- 0.05/num.obs
		Bonf <- rep(FALSE, num.obs)
		Bonf[p.valC<Bonf.CO] <- TRUE
	
		FDR.01 <- rep(FALSE, num.obs)
		FDR.01[bs$qval<0.01] <-  TRUE

		FDR.05 <- rep(FALSE, num.obs) 
		FDR.05[bs$qval<0.05] <-  TRUE
	
		out <- cbind(bs, tail=Tail, p.val.tail=p.valC, Bonf, FDR.01, FDR.05)
		write.table(out, file=paste(filepath, ".Cpval", sep=""))
		return(out)
}
