#dependencies and all source files
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 

#make sure you have latest package of hierfstat: version (0.04-10) earlier versions give errors
#http://cran.r-project.org/web/packages/hierfstat/index.html

if (is.installed("hierfstat")==FALSE){
  install.packages("hierfstat",type="source", dependencies=TRUE)
}

if (is.installed("SuppDists")==FALSE){
  install.packages("SuppDists",type="source", dependencies=TRUE)
}

if (is.installed("ggplot2")==FALSE){
  install.packages("ggplot2",type="source", dependencies=TRUE)
}

if (is.installed("lattice")==FALSE){
  install.packages("lattice",type="source", dependencies=TRUE)
}

library(hierfstat)
library(SuppDists)
library(ggplot2)
library(lattice)

## install Storey's q-value package
source("http://bioconductor.org/biocLite.R")
biocLite("qvalue")
library(qvalue)