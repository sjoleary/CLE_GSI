
getAlleles <- function(column, diploid, ndig){
### This function gets alleles from a vector of diploid genotypes in FSTAT format
### diploid is a logical TRUE or FALSE
### ndig is the number of digits in the coding for alleles
### returns a vector of alleles  
  
  if (diploid==TRUE){
			if (ndig==1){
				b <- cbind(substr(as.character(column), start=1,stop=1), substr(as.character(column), start=2,stop=2))
				}
			if (ndig==2){
				b <- cbind(substr(as.character(column), 1,2), substr(as.character(column), 3,4))
				}	
			 if (ndig==3){
				b <- cbind(substr(as.character(column), 1,3), substr(as.character(column), 4,6))
				}
		}	
		if (diploid==FALSE){	
				b <- column												 
		}	
	return(b)
}

getAllCounts <- function(pops, column, diploid, ndig){
### This function returns allele counts in each population
### pops is a vector assigning each individual to each population
### column is a vector of diploid genotypes
### diploid is a logical
### ndig is the number of digits used to code the genotypes
	a<- as.vector(getAlleles(column, diploid, ndig))
	if (diploid==TRUE){b<-tapply(a, list( c(pops,pops), a), length)}else{
	b<-tapply(a, list(pops,a), length)}
	b[which(is.na(b))] <- 0
	return(b)
	}