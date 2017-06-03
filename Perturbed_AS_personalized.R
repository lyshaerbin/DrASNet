####This script is used to identify the personalized
####alternative splicing perturbations in cancer
###Input:
###PSImat:The filename that store the PSI matrix for cancer and normal samples
###PSIsam:The sample file for the tumor and normal
Personalized_AS=function(PSImat,PSIsam){
  cancerAS <- read.csv(paste(PSImat,".txt",sep = ""),stringsAsFactors = FALSE, header = FALSE,sep = "\t")
  C_sample<-read.csv(paste(PSIsam,".txt",sep = ""),stringsAsFactors = FALSE, header = FALSE,sep = "\t")
  nor_Sample<-grep("_Norm",C_sample$V1)
  can_Sample<-setdiff(c(1:dim(C_sample)[1]),nor_Sample)
  N<-dim(cancerAS)[2]
  Cancerexp<-cancerAS[,(11:N)]
  CancerOUT<-c()
  for (j in 1:dim(cancerAS)[1]) {
    print(j)
    outlier <- boxplot.stats(as.numeric(Cancerexp[j,can_Sample]))$out
    outlier<-unique(outlier)
    temps<-c()
    KW<-c()
    if (length(outlier)>0){
      for (k in 1:length(outlier)) {
        xx<-which(Cancerexp[j,can_Sample]==outlier[k])
        for (mm in 1:length(xx)) {
          temps<-rbind(temps,C_sample$V1[can_Sample[xx[mm]]])
        }
      }
      KW<-cbind(temps,cancerAS[j,1:3])
    }
    CancerOUT<-rbind(CancerOUT,KW)
  }
  colnames(CancerOUT)=c("SampleID","ASgene","ASID","AStype")
  return(CancerOUT)
}



