#########This script prioritized the trans mutation-AS pairs
Trans_pri=function(TRANS_mut){

  mut_AS<-unique(TRANS_mut[,2:3])

  xgene<-unique(mut_AS[,1])

  Mutlist<-c()
  while(dim(mut_AS)[1]>0){
    ###compute the degree of each mutated gene
    gene<-unique(mut_AS[,1])
    degree<-c()
    for (i in 1:length(gene)) {
      degree<-rbind(degree,length(which(mut_AS[,1]==gene[i])))
    }
    maxdeg<-max(degree)
    temgene<-which(degree==maxdeg) ###tha max degree
    Mutlist<-rbind(Mutlist,gene[temgene[1]]) ###add the mut gene to the list
    xx<-which(mut_AS[,1]==gene[temgene[1]]) ###find the AS explained by mut gene
    exAS<-mut_AS[xx,2]
    mut_AS<-mut_AS[-xx,]
    kx<-c()
    for (j in 1:length(exAS)) {
      nn<-which(mut_AS[,2]==exAS[j])
      print(nn)
      if (length(nn)>0){
        kx<-rbind(kx,nn)
        #print(kx)
      }

    }
    if(length(kx)>0){
      mut_AS<-mut_AS[-kx,]
    }
    #print(dim(mut_AS))
  }
  Mutlist<-unique(Mutlist)
  return(Mutlist)
}
###############example

