change.extension <- function(infilepath, newExt){
  s1 <- unlist(strsplit(x=infilepath, split="\\/")) #get the filename from the path
  filename <- s1[length(s1)]
  s2 <- unlist(strsplit(x=filename, split="\\."))#separate extension from filename
  outfile = paste(paste(s2[-length(s2)], collapse=".", sep=""), newExt, sep="")
  outfilepath = paste(paste(s1[-length(s1)],collapse="/"), outfile, sep="/")
  return(outfilepath)
}