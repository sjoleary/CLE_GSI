returnTimestamp <- function(){
  round(unclass(Sys.time()))
}

getRunName <- function(){
  return(paste(round(unclass(Sys.time())), sample(1001:9999,1), sep="."))
}
