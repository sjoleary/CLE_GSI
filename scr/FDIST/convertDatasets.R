
ConvertGenepopToFSTAT <- function(GenepopFile, NumLoci, HighestLociValue, 
                                  NumDigitsLoci, OutFilePath){
### Function for converting a datafile from strict Genepop format to FSTAT format
### First line of FSTAT file contains 4 numbers: 
  # np: the number of sampled populations, 
  # NumLoci: the number of loci, 
  # HighestLociValue: the highest number used to label an allele
  # NumDigitsLoci: 
      # a number for the 1 if the code for alleles is a one digit number (1-9), 
      # a 2 if code for alleles is a 2 digit number (01-99) or a 
      # a 3 if code for alleles is a 3 digit number (001-999). 
  # These 4 numbers need to be separated by any number of spaces.

### The first line is immediately followed by nl lines, 
    # each containing the name of a locus, 
    # in the order they will appear in the rest of the file.
	
	file <- read.table(GenepopFile, header=FALSE, sep="/")
	lineIDs <- which(file[,1]=="pop")
	lineIDs <- c(lineIDs, dim(file)[1]+1)
	locinames <- file[2:(lineIDs[1]-1),]
	if (NumLoci!=length(locinames)){print("Number of loci is not equal to number of loci names in file")}
	popID <- vector()
	for (i in 2:length(lineIDs)){
			popID <- c(popID, rep(i-1, (lineIDs[i] - lineIDs[i-1] -1)))
	}
	length(popID)
	FSTATdata <- file[lineIDs[1]:dim(file)[1],]
	FSTATdata <- FSTATdata[-(which(FSTATdata=="pop"))]
	length(FSTATdata)
	FSTATdata2<- substring(FSTATdata, 5, nchar(as.character(FSTATdata)))
	FSTATdata3 <- paste(popID, FSTATdata2, sep="")
	
	NumSamples <- length(FSTATdata2)
	write(paste(NumSamples, NumLoci, HighestLociValue, NumDigitsLoci, sep=" "), OutFilePath, ncol=1)
	write(as.character(locinames), OutFilePath, ncol=1, append=TRUE)
	write(FSTATdata3, OutFilePath, ncol=1, append=TRUE)
}