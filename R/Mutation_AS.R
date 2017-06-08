####This script is used to identify the personalized Mutation-
####alternative splicing perturbations in cancer
###Input:
###Cancer: The tumor type
###SampleAS: Personalized AS events
###mutation: Personalized mutations
###network: Protein-protein interaction network
Mutation_AS=function(SampleAS,mutation,network){
  AS<-SampleAS
   # the unique samples in each cancer
  U_sample_mut<-unique(mutation$V1)
  U_sample_AS<-unique(AS[,1])
  U_sample<-intersect(U_sample_AS,U_sample_mut)

  cIS_MUT<-c()
  TRANS_MUT<-c()

  for (i in 1:length(U_sample)) {
    x<-which(mutation$V1==U_sample[i])
    Smut<-mutation$V3[x]
    Smut<-unique(Smut)
    y<-which(AS[,1]==U_sample[i])
    SAS<-AS[y,2]
    SAS<-unique(SAS)
    for (j in 1:length(Smut)) {
      for (k in 1:length(SAS)) {
        print(c(i,j,k))
        gene1<-Smut[j]
        gene2<-SAS[k]
        if(identical(gene1,gene2)==T){
          cIS_MUT<-rbind(cIS_MUT,cbind(U_sample[i],gene1,gene2)) ###cis regulation
        } else{
          kk<-which(network$V1==gene1&network$V2==gene2)
          if (length(kk)>0){
            TRANS_MUT<-rbind(TRANS_MUT,cbind(U_sample[i],gene1,gene2)) ###trans regulation
          }
        }
      }
    }
  }
  Mut.AS=list(Cis=cIS_MUT,trans=TRANS_MUT)
  return(Mut.AS)
}

#########################
##examples
# mutation<-read.csv("KICH_sample_mut_code.txt",stringsAsFactors = FALSE, header = FALSE,sep = "\t")
# AS<-read.csv("KICHoutlier.txt",stringsAsFactors = FALSE, header = FALSE,sep = "\t")
# network<-read.csv("network.txt",stringsAsFactors = FALSE, header = FALSE,sep = "\t")
# MAS=Mutation_AS(AS,mutation,network)
