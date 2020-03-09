


getHe.dataset <- function(data1, diploid, ndig=NULL){
########################   
##### Get He for each locus from a dataset in FSTAT format
########################  
	numloci <- ncol(data1)-1
	He <- rep(NA, numloci)
  lociMat <- data1[,-1]
  He <- apply(lociMat,2, getHe.locus, diploid, ndig)
	return(He)
}

getHe.locus <- function(vector, diploid, ndig=NULL){
########################   
##### Get He for a single locus column from a dataset
########################    
  b <- getAlleles(vector, diploid=diploid, ndig=ndig)														 																
	b2 <- unlist(b)
		
  ### remove NAs if they exist
  	whichNA <- which(is.na(b2)==TRUE)
  	if (length(whichNA)>0){
      b2 <- b2[-whichNA]
  	}	
		
  a <- tapply(b2, b2, length)
	He <- 1 - sum((a/sum(a))^2)
  return(He)
}