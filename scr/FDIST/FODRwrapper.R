### FODR wrapper
### 1) Set FST calculation and upload data
### 2) run IM sims
### 3) Upload sims and evaluate p-values
### 4) Make FODR plots

### 1) Set FST calculation
  FST.set("WC84")

### 2) Upload data in FSTAT format
  infile <- "data/FSTAT_1R_Neut9000_75_20_Set1.dat"
  data1 <- read.fstat(infile)
  ### recommend that data is in actual FSTAT format, because it specifies the 
  ### number of alleles used for coding, compatible with the hierfstat package

### 3) calculate FST for each locus and overall FST
  data1.FSTs <- FSTcalc.dataset(data1, diploid=FALSE, ndig=1)
  str(data1.FSTs)
  
### 4) Run IM simulations 
  ndemes <- 10 ### specify the number of demes. 
    ### What you choose for ndemes may largely determine the # outliers, because
    ### it determines the variance in the FST distribution. Especially when ndemes < 10
  FST.overall <- data1.FSTs$FST.overall
    ### overall FST from the data
  sample.size <- get.corrected.mean.sample.size(data1) 
    ### Specify the sample sizes for each dems in the IM simulations
    ### In this case the sample sizes are equal from each deme, 
    ### and based on a weighted mean of the sample sizes in the data
    ### It is also possible to input a vector here, but the length MUST be equal to ndemes
  
  outfilepath.sims <- change.extension(infile, "All.sims")
  number.reps <- 100000
    #expect ~5min per 1000 reps. of 20-deme IM.  so 8 hours for 100,000 reps
  
  IM.sims.replicate(ndemes, FST.overall, sample.size, number.reps, outfilepath.sims, 
                  diploid=TRUE, gens=NULL, append.outfile=FALSE, progress=TRUE)
  
### 5) Input the simulation file
  FST.sim <- read.table(outfilepath, header=TRUE)
  
### 6) Evaluate p-values for each locus in the dataset and write to file
  FST.emp <- data1.FSTs$FST.emp
  outfilepath.pvals <- change.extension(infile, ".pval")
  All.Pvals <- write.all.PvalueInfo.dataset(FST.emp, FST.sim, outfilepath.pvals)
  
