#######################################################################  	
getIM.Params.Haploid <- function(totNumDemes, obsFST){
    #this is for a finite number of demes
    N=1000
    m=(1.0/obsFST-1.0)*(totNumDemes-1.0)/(2.0*N*totNumDemes)
    return(list(ndemes=totNumDemes, m=m, N=N))
}
#######################################################################  
getIM.Params.Diploid <- function(totNumDemes, obsFST){
    #this is for a finite number of demes
    #assume N
    N=1000
    m=(1.0/obsFST-1.0)*(totNumDemes-1.0)/(4.0*N*totNumDemes)    
    return(list(ndemes=totNumDemes, m=m, N=N))
}

#######################################################################  
loadCFunctions.withDelete <- function(sourceWD){		
		setwd(sourceWD)	
		if(file.exists("src/SimIM.o")){system(paste("rm src/SimIM.o"))}
		if(file.exists("src/SimIM.so")){system(paste("rm src/SimIM.so"))}
		system("R CMD SHLIB src/SimIM.c")
		dyn.load("src/SimIM.so")
}
loadCFunctions <- function(sourceWD){
	system("R CMD SHLIB src/SimIM.c")
	dyn.load("src/SimIM.so")
}

#######################################################################  
simIM.neut.return.p <- function(ndemes, N, m, u=NULL, p.start=runif(1), gens)  {	
########################
### This function runs one replicate simulation of the island model and 
  ### returns the allele frequency in each deme
### ndemes is the number of demes in the IM
### N is the population size per deme
### m is the migration rate between demes
### default mutation rate (u) is set to m/1000 because we found this resulted in
  ### more or less stable dynamics for a wide range of N and m
### p.start is the starting allele frequency in all demes
  ### randomly chosen between 0 and 1 
  ### (this results in a wide spread of He over several replicates)
### gens is the number of generations

		if(length(u)==0){
			u <-  m/100
		}
	
		p.freq <- rep(p.start, ndemes)

		out.p <- .C("SimulateIM", as.double(p.freq), as.integer(length(p.freq)), 
                as.double(u), as.double(m), as.integer(N), as.integer(gens)) 
		p.freq.out <- out.p[[1]]
		
		return(list(p.start=p.start, p.freq.out=p.freq.out, u=u)) #3 columns		
}

#######################################################################  
simIM.selection.return.p <- function(ndemes, N, m, u=NULL, 
                                     p.start=runif(1), gens, S.vect)	{	
########################
### This function runs one replicate simulation of the island model with selection
### ndemes is the number of demes in the IM
### N is the population size per deme
### m is the migration rate between demes
### default mutation rate (u) is set to m/1000 because we found this resulted in
  ### more or less stable dynamics for a wide range of N and m
### p.start is the starting allele frequency in all demes (single value)
  ### randomly chosen between 0 and 1 
  ### (this results in a wide spread of He over several replicates)
### gens is the number of generations
### S.vect is the vector of selection coefficients in each deme

	if(length(u)==0){
		u <-  m/1000
	}
	
	p.freq <- rep(p.start, ndemes)
	
	out.p <- .C("SimulateIM.Selection", as.double(p.freq), 
              as.integer(length(p.freq)), as.double(u), as.double(m), 
              as.integer(N), as.integer(gens), as.double(S.vect)) 

	p.freq.out <- out.p[[1]]
	return(list(p.start=p.start, p.freq.out=p.freq.out, u=u))
}
#######################################################################  
get.selVECT <- function(s.low,s.high, envi.vect){
########################
### This function gets selection coefficients from an environmental vector
### Assumes a linear relationship between environment and selection
### s.low and s.high are approximately the lowest and highest coefficients on the landscape
  if (s.high==0 & s.low==0){s.VECT=0}else{
		if(s.high==-s.low){
			s.slope <- (s.high - s.low)/(4*sd(envi.vect))
        #slope is determined by 4 SD (ie ~ 95% CI)
        #intercept is 0 when s.low=s.high
			s.vect <- s.slope*envi.vect
		}else{
			s.slope <- (s.high - s.low)/(4*sd(envi.vect))
			s.int <- s.high-2*sd(envi.vect)*s.slope 
			s.vect  <- s.slope*envi.vect + s.int
		} #end nested if
	}		#end if
	return(s.vect)
}

#######################################################################  
sample.pfreq <- function(p.freq.out, sample.sizes){
### This function takes a vector of true allele frequencies and returns the 
    ### counts of alleles sampled in each population (alleles in columns, pops in rows)
### sample.sizes is a vector of sample sizes for each population 
    ### can be a single value
    ### or a vector (does not have to be equal for each deme)
### For haploid sampling
    ndemes <- length(p.freq.out)
    if (length(sample.sizes)>1  & length(sample.sizes)!=ndemes){
      print("Error in sample.pfreq: sample size vector does not equal number of demes in IM")
      break;
    }
		samples <- rbinom(ndemes, sample.sizes, p.freq.out)
	return(matrix(c(samples, sample.sizes-samples), ncol=2))	
}

#######################################################################  
do.IM.sim.get.FST <- function(ndemes, N, m, sample.sizes, 
                              diploid=TRUE, gens=NULL){
########################   
### This is a wrapper for one locus, to do an island model simulation and return FST

  if (length(gens)==0){gens <- get.num.gens(m, N)}
    # If number gens not specified, run until half-life

  if (length(sample.sizes)>1 & length(sample.sizes)!=ndemes){
      print("Error in sample.pfreq: sample size vector does not equal number of demes in IM")
    }
  
  p.freqs <- simIM.neut.return.p(ndemes, N, m, gens=gens)
  alleleCounts <- sample.pfreq(p.freqs$p.freq.out, sample.sizes)
  return(FSTcalc.sims(alleleCounts))

}

#######################################################################  
IM.sims.replicate <- function(ndemes, FST.overall, sample.sizes, number.reps, 
                              outfilepath, diploid=TRUE, gens=NULL, 
                              append.outfile=FALSE, progress=TRUE){
########################   
### Do replicate simulations of the island model and writes to file
### Also writes a metadata file for these simulations
### every 1000 runs it appends the list of FSTs to "outfilepath"
### ndemes is the number of demes in the IM
### FST.overall is the mean FST for the IM simulations
### sample.sizes is a vector of sample sizes for each population 
    ### can be a single value
    ### or a vector (does not have to be equal for each deme)
### number.reps is the number of replicate loci to run with these parameters

  if(file.exists(outfilepath) & append.outfile==FALSE){
    writeLines("Your outfile already exists.  If you would like to append these
               simulations to this outfile, specify append.outfile=TRUE")
    break;
  }

  if(file.exists(outfilepath)==FALSE & append.outfile==TRUE){
    writeLines("You specified that you want to append to an outfile, but that 
                outfile does not exist.")
    break;
  }

  if (diploid==TRUE){
		IMParams <- getIM.Params.Diploid(ndemes, FST.overall)
  }else{
    IMParams <- getIM.Params.Haploid(ndemes, FST.overall)
  }

  ## Write information to a metafile
  metafileName <- change.extension(outfilepath, ".meta")
  #print(metafileName)
  write(paste("*************************************"), file=metafileName, append=TRUE)
  write(paste("Simulations started on ", Sys.time()), file=metafileName, append=TRUE)
  write(paste("Simulations writing to file: ", getwd(), "/", outfilepath, sep=""), 
        file=metafileName, append=TRUE)

  if (length(gens)==0){
    write("Number of generations chosen automatically based on half-life", 
          file=metafileName, append=TRUE)  
    gens <- get.num.gens(IMParams$m,  IMParams$N)
  }
  
  write(paste("Number of generations:", gens), file=metafileName, append=TRUE)
	
  if (progress==TRUE){
    writeLines(paste("Starting the simulations... 
Will append to simulation file with the name: ", outfilepath))
  }

#  system.time(
    # Expect about 5 minutes for every 1000 samples
    # 25 minutes for every 5000 samples
 #   t(replicate(10, do.IM.sim.get.FST(ndemes=20, N=1000, m=0.005, sample.sizes=30)))
#  )
tot<-0
x <- 0
per.batch=1000
TOT <- 0

    while (TOT < number.reps){
      if(TOT<(number.reps-per.batch)){
      out <- t(replicate(per.batch, 
                         do.IM.sim.get.FST(ndemes=IMParams$ndemes, 
                                           N=IMParams$N, 
                                           m=IMParams$m, sample.sizes=sample.sizes)))
      }else{
        out <- t(replicate((number.reps-TOT), 
                            do.IM.sim.get.FST(ndemes=IMParams$ndemes, 
                                              N=IMParams$N, m=IMParams$m, 
                                              sample.sizes=sample.sizes)))
      }
      out2 <- as.data.frame(matrix(as.numeric(out), ncol=ncol(out)))
        #out is a list.  Have to transform to numeric and dataframe
      colnames(out2)<- colnames(out)
      out3 <- out2[complete.cases(out2),]
        #remove loci that drifted to fixation
    
      if (file.exists(outfilepath)==FALSE){
        write(colnames(out3), file=outfilepath, ncolumns=ncol(out3))
      }
    
      x <- x + nrow(out2)
      tot<- tot + nrow(out3)
        #keep track of total number of loci simulated
      if (progress==TRUE){
        writeLines(paste("Kept", tot, "of", x, "runs",sep=" "))
      }
      write(paste("Kept", tot, "of", x, "runs",sep=" "), file=metafileName, append=TRUE)
    
      write.table(out3, file=outfilepath, row.names=FALSE, col.names=FALSE, append=TRUE)
      TOT <- tot
    }#end while loop
  write(paste("----------------------------------"), file=metafileName, append=TRUE)
  write(paste("Simulations ended on ", Sys.time()), file=metafileName, append=TRUE)
  write(paste("Total number simulations:", TOT ), file=metafileName, append=TRUE)
  write(paste("Overall FST parameter:", FST.overall), file=metafileName, append=TRUE)
  write(paste("Island Model parameter ndemes:", IMParams$ndemes), 
        file=metafileName, append=TRUE)
  write(paste("Island Model parameter N:", IMParams$N), file=metafileName, append=TRUE)
  write(paste("Island Model parameter m:", IMParams$m), file=metafileName, append=TRUE)
  write(paste("Number of replicates requested:", number.reps), 
        file=metafileName, append=TRUE)
  write(paste("Vector of sample sizes from each deme:"), file=metafileName, append=TRUE)
  write(sample.sizes, file=metafileName, append=TRUE, ncol=length(sample.sizes))
  write(paste("***********************************"), file=metafileName, append=TRUE)

}
