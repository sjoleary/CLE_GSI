#library(klaR)
library(MASS)
library(MVN)

# Assess the accuracy of the prediction
# percent correct for each category of G

summarize.class<-function(original, classify) {
  class.table<<-table(original, classify)
  numb<<-rowSums(class.table)
  prop<<-round(class.table/numb,4)
  list(class.table = class.table, prop = prop)
}

BoxCox <- function(data, type = c("optimal", "rounded")){
  
  powerTransformation = summary(powerTransform(data))$result
  
  if(type == "optimal"){
    
    lambda = powerTransformation[,1]
    
  }
  
  if(type == "rounded"){
    
    lambda = powerTransformation[,2]
    
  }
  
  for(i in 1:length(lambda)){
    
    if(lambda[[i]] == 0){
      
      data[i] = log(data[i])
      
    }else{
      
      data[i] = data[i]^lambda[[i]]
      
    }
  }
  
  result = list(data, lambda)
  
  return(result)
}

boot.qda<-function(data, indices){
  tryCatch(
    {
      data<-data[indices,]
      shmerl<<-data
      x<<-qda(as.formula(as.character(result[1,1])), data=data ,CV=TRUE,lambda=0,gamma=0)
      y<<-summarize.class(data[,names(data)==groupvar], x$class)
      return(c(diag(y$prop),mean(diag(y$prop))))
    },
    error=function(cond){
      foo<<-0
      return(foo)
    }
  )
}

successest<-function(x){
  tryCatch(
    {
      foo<<-x
      test<<-as.formula(foo[1])  
      fit <<- qda(test, data=dataname ,CV=TRUE,lambda=0,gamma=0)
      class.succ<<-summarize.class(dataname[,names(dataname)==groupvar], fit$class)
      # Calculate Klecka's Tau
      n.corr<<-sum(diag(class.succ$class.table))
      nipi<<- apply(class.succ$class.table,MARGIN=1,FUN=sum)*(1/ncol(class.succ$class.table))
      numerator<<-n.corr-sum(nipi)  
      N<<-sum(class.succ$class.table)
      denominator<<-N-sum(nipi)
      tau<<-numerator/denominator
      
      # Calculate overall classification success (% of all samples correctly classified)
      Over.Succ<<-sum(diag(class.succ$class.table))/sum(class.succ$class.table)
      foo[2]<<-as.numeric(Over.Succ)
      foo[3]<<-as.numeric(tau)
      return(foo)
    },
    error=function(cond){
      foo[2]<<-0
      return(foo)
    }
    
  )
}

# Select grouping here by un-commenting: Estuary","ThreeRegion","TwoRegion"
groupvar<-"Estuary"
#groupvar<-"ThreeRegion"
#groupvar<-"TwoRegion"

# Select year targeted for analysis by un-commenting below
year<-2013
#year<-2014
#year<-2015
#year<-2016

#Set working directory to location of chemistry data .csv prior to reading in.
dat<-read.csv("Cleucas_chemistry_data.csv")
dat$Estuary<-factor(dat$Estuary,levels(dat$Estuary)[c(4,2,3,1)])
dat$ThreeRegion<-factor(dat$ThreeRegion,levels(dat$ThreeRegion)[c(2,1,3)])

dat<-dat[,-14]
dat$Sampled.Year<-as.factor(dat$Sampled.Year)
dat$Year.Collected<-as.factor(dat$Year.Collected)

dat<-dat[dat$Sampled.Year==year,]
dat$Sampled.Year<-droplevels(dat$Sampled.Year)

#make values non-negative
dat$d13C<-dat$d13C+15
dat$d18O<-dat$d18O+5

row.names(dat)<-seq(1,nrow(dat),1)


MVnorm1<-mvn(dat[,c(10:18)],mvnTest="hz",bc=TRUE,multivariatePlot = "qq",multivariateOutlierMethod = "adj",showNewData=TRUE,showOutliers = TRUE)
MVnorm1$multivariateNormality

dat3<-dat

##########################################
dat2<-dat[row.names(MVnorm1$newData),]

newdat<-MVnorm1$newData
newdat$mergecols<-as.numeric(row.names(newdat))

metadat<-dat2[,c(1:9)]
dat3<-cbind(metadat,newdat)

###########################################
dat3$SharkID<-droplevels(dat3$SharkID)

vars<-c("d13C","d18O","Li","Mn","Mg","Sr","Cu","Ba")
combinations<-vector("list")

for(i in 0:(length(vars)-1)){
  a<-combn(vars,length(vars)-i)
  b<-apply(a,2,paste,collapse="+")
  #var.df[i+1,1]
  combinations[[i+1]]<-b
}

frmlas<-unlist(combinations)

formulas<-paste(groupvar,"~",frmlas,sep=" ")
successcombos<-data.frame(Formulas=formulas,Success=rep("X",length(formulas)),Tau=rep("X",length(formulas)))
successcombos[,1]<-as.character(successcombos[,1])

dataname<-dat3

result<-as.data.frame(t(apply(successcombos,1,FUN=successest)))
result[,2]<-as.numeric(as.character(result[,2]))
result<-result[order(-result$Success),]
names(result)[3]<-"Tau"
View(result)[1]


#get QDA result from combination of variables producing highest overall classification success
bestfit<-qda(as.formula(as.character(result[1,1])), data=dataname ,CV=TRUE,lambda=0,gamma=0)
