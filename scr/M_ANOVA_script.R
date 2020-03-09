# Set working directory to location of vert chemistry data .csv before reading in.
dat<-read.csv("Cleucas_chemistry_data.csv")
dat$Estuary<-factor(dat$Estuary,levels(dat$Estuary)[c(4,2,3,1)])
dat<-dat[,-14]

dat13<-dat[dat$Sampled.Year==2013,]
dat13<-dat13[dat13$Estuary!="Matagorda",]
dat13$Estuary<-droplevels(dat13$Estuary)
dat14<-dat[dat$Sampled.Year==2014,]
dat14$Estuary<-droplevels(dat14$Estuary)
dat15<-dat[dat$Sampled.Year==2015,]
dat15$Estuary<-droplevels(dat15$Estuary)
dat16<-dat[dat$Sampled.Year==2016,]
dat16$Estuary<-droplevels(dat16$Estuary)

yearaov<-function(X,DATA){
  form<-as.formula(paste(TRACER,"~Estuary"))
  summary(aov(form,data=DATA))
}

results<-data.frame(matrix(ncol=4,nrow=9))
names(results)<-c("Tracer","DF","Fval","p")
j=0
for (i in 10:18){
  j<-j+1
  TRACER<-names(dat13)[i]
  form<-as.formula(paste(TRACER,"~Estuary"))
  a<-summary(aov(form,data=dat13))
  results[j,1]<-TRACER
  results[j,2]<-paste(a[[1]][[1]][1],a[[1]][[1]][2],sep=",")
  results[j,3]<-a[[1]][[4]][1]
  results[j,4]<-a[[1]][[5]][1]
  }


tracers<-as.matrix(cbind(dat[,10:18]))

# Run two-way MANOVA with interaction
summary(manova(tracers~Estuary*Year.Collected,data=dat))


# Run one-way ANOVAs with interactions
summary(aov(Li~Estuary*Year.Collected,data=dat))
summary(aov(d13C~Estuary*Year.Collected,data=dat))
summary(aov(d18O~Estuary*Year.Collected,data=dat))
summary(aov(Mg~Estuary*Year.Collected,data=dat))
summary(aov(Mn~Estuary*Year.Collected,data=dat))
summary(aov(Cu~Estuary*Year.Collected,data=dat))
summary(aov(Zn~Estuary*Year.Collected,data=dat))
summary(aov(Sr~Estuary*Year.Collected,data=dat))
summary(aov(Ba~Estuary*Year.Collected,data=dat))

