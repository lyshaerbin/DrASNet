##This script identify the perturbed alternative splicing
##events in cancer
###Input:
###PSImat:The filename that store the PSI matrix for cancer and normal samples
###PSIsam:The sample file for the tumor and normal
###FDRthr: The threshold of FDR
###FCthr: The threshold of fold change
PerturbeAS=function(PSImat,PSIsam,FDRthr,FCthr){
  cancerAS <- read.csv(paste(PSImat,".txt",sep = ""),stringsAsFactors = FALSE, header = FALSE,sep = "\t")
  C_sample<-read.csv(paste(PSIsam,".txt",sep = ""),stringsAsFactors = FALSE, header = FALSE,sep = "\t")
  nor_Sample<-grep("_Norm",C_sample$V1)
  can_Sample<-setdiff(c(1:dim(C_sample)[1]),nor_Sample)
  N<-dim(cancerAS)[2]
  if (length(nor_Sample)>=5){
    FC<-c()
    P<-c()
    Cancerexp<-cancerAS[,(11:N)]
    print("begin to calculated the AS events:")
    for (j in 1:dim(cancerAS)[1]) {
      print(j)
      mec<-mean(as.numeric(Cancerexp[j,can_Sample]))/mean(as.numeric(Cancerexp[j,nor_Sample]))
      ptemp<-wilcox.test(as.numeric(Cancerexp[j,can_Sample]),as.numeric(Cancerexp[j,nor_Sample]))$p.value
      FC<-rbind(FC,mec)
      P<-rbind(P,ptemp)
    }
    result<-cbind(cancerAS[,1:10],FC[,1],P[,1])
    fdr<-p.adjust(as.numeric(result$`P[, 1]`),method = "BH")
    x<-which(as.numeric(result$`FC[, 1]`)>FCthr&fdr<FDRthr)
    y<-which(as.numeric(result$`FC[, 1]`)<(1/FCthr)&fdr<FDRthr)
    z<-union(x,y)
    KK<-cbind(result[z,],fdr[z])
    colnames(KK)=c("ASgene","ASid","AS.type","Exons","From.exon",
                   "Toexon","Novel","pct_with_values","psi_range",
                   "std_psi","Fold.change","P.value","FDR")
    return(KK)
  }
}

#############Using CHOL cancer ans an example
PSImat="CHOL_mean"
PSIsam="CHOL_sample"
KW=PerturbeAS(PSImat,PSIsam,0.01,2)
