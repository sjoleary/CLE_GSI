getCI <- function(vector, percentile){
  a<-sort(vector)
  index2=round(length(a)*percentile,0)
  index1=round(length(a)*(1-percentile),0)
  return(c(a[index1], a[index2]))
}
get.CI.He <- function(xseq, percentile, FSTs.sim) {
  CI.He<-  matrix(NA, length(xseq),2)
  for (i in 1:(length(xseq)-1)){
     indexes <- which(xseq[i]<FSTs.sim$He & xseq[i+1]>=FSTs.sim$He)
     if(length(indexes>0)){
     CI.He[i,]<-getCI(FSTs.sim$FST[indexes], percentile)
     }
  }
  return(CI.He)
}

makeCIplot <- function(percentile,c, FSTs.sim){
  CI.He<-get.CI.He(xseq, percentile, FSTs.sim)
  points(xseq, CI.He[,1], type="l", col=c, lwd=2)
  points(xseq, CI.He[,2], type="l", col=c, lwd=2)
}



make.He.FST.plot <- function(FSTs.sim, All.Pvals, truepos, FODR.my.plot=FALSE){
  maxHe <- max(c(FSTs.sim$He, All.Pvals$He))
  minHe <- min(c(FSTs.sim$He, All.Pvals$He))
  maxFST <- max(c(FSTs.sim$FST, All.Pvals$FST))
  minFST <- min(c(FSTs.sim$FST, All.Pvals$FST))
  
  He.cell <- 0.01
  
  xseq<-seq(0,maxHe+He.cell,He.cell)
  yseq<-seq(minFST-FST.cell ,maxFST+FST.cell ,FST.cell )
  
  plot(NULL, NULL, xlim=c(0,maxHe), ylim=c(min(yseq), max(yseq)), 
       bty="l", las=1, 
       ylab=expression(F[ST]), xlab=expression(H[e]))
  
  points(All.Pvals$He, All.Pvals$FST, col=rgb(.5,.5,.5,0.1), pch=19, cex=0.8)
  points(All.Pvals$He[truepos], All.Pvals$FST[truepos], col="black", pch=19, cex=0.8)
  
  if(FODR.my.plot==TRUE){
    makeCIplot(0.95,"red", FSTs.sim)
    makeCIplot(0.96,"orange", FSTs.sim)
    makeCIplot(0.97,"yellow", FSTs.sim)
    makeCIplot(0.98,"green", FSTs.sim)
    makeCIplot(0.99,"blue", FSTs.sim)
    makeCIplot(0.999,"purple", FSTs.sim)
  }else{
    makeCIplot(0.95, "blue", FSTs.sim)
    makeCIplot(0.999,"blue", FSTs.sim)
  }
}






# br<-c(-0.5:10, seq(10,200, by=10), seq(100,3000, by=100))
# image(xseq[1:51]+He.cell/2,yseq[1:51]+FST.cell/2,
#       (z[1:51,1:51]), zlim=c(1,max(z)),
#       breaks = br,
#       col=c(rev(grey.colors((length(br)-1),  start = 0.2, end = 0.5))),
#       xlim=c(0,1),ylim=c(minFST, maxFST)      
#       )
# 
# plot(NULL, NULL, xlim=c(0,0.55), ylim=c(0,0.1))
# 
# z<- matrix(NA, length(xseq)-1, length(yseq)-1)
# for (i in 1:(length(xseq)-1)){
#   for (j in 1:(length(yseq)-1)){
#     #would be faster to get counts on histogram every step
#     index <- xseq[i]<FSTs.sim$He & xseq[i+1]>=FSTs.sim$He & yseq[j]<FSTs.sim$FST & yseq[j+1]>=FSTs.sim$FST
#     z[i,j] <- sum(index)
#     
#   }
# }

#getCount <- function(i,j){
#  index <- xseq[i]<FSTs.sim$He & xseq[i+1]>FSTs.sim$He & yseq[j]<FSTs.sim$FST & yseq[j+1]>FSTs.sim$FST
#}
#a<-mapply(getCount,xseq,yseq)



