get.sample.size.each.pop <- function(data1){
  #data1 is in FSTAT format
  #assume 1st column is population of that individual
  as.numeric(unlist(tapply(data1[,1], data1[,1], length)))
}

get.mean.sample.size <- function(data1){
  #data1 is in FSTAT format
  #assume 1st column is population of that individual
  mean(as.numeric(unlist(tapply(data1[,1], data1[,1], length))))
}

get.corrected.mean.sample.size <- function(data1){
  sample.sizes <- as.numeric(unlist(tapply(data1[,1], data1[,1], length)))
  n.pops <- length(sample.sizes)
  n.ave <- mean(sample.sizes)
  n.c <- (n.pops*n.ave - sum(sample.sizes^2)/(n.pops*n.ave))/(n.pops-1) 
  return(n.c)
}