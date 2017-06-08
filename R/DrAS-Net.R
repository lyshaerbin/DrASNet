#####DrAS-Net main funciton
#################using the KICH as an example
####identify the Personalized Alternative splicing events in cancer
PSImat="CHOL_mean"
PSIsam="CHOL_sample"
source("Perturbed_AS_personalized.R")
#####this process might take some time if there are many AS events
Person.AS=Personalized_AS(PSImat,PSIsam)
write.table(Person.AS,"Person_AS_event.txt",row.names=F,col.names=T,sep="\t",quote=F)
################
##Identify the mutation-AS pairs in all samples
###read the mutation file
mutation<-read.csv("CHOL_sample_mut_code.txt",stringsAsFactors = FALSE, header = FALSE,sep = "\t")
###read the network file
network<-read.csv("network.txt",stringsAsFactors = FALSE, header = FALSE,sep = "\t")
source("Mutation_AS.R")
#####identify the mutation-AS pair
#####this process might take some time if there are many AS events and samples
Per.M.AS=Mutation_AS(Person.AS,mutation,network)
######TRANS mutation-AS pairs
Trans_mut=Per.M.AS$trans
write.table(Trans_mut,"Trans_mut_AS.txt",row.names=F,col.names=F,sep="\t",quote=F)

###trans gene priorization
source("Trans_pri.R")
Trans=Trans_pri(Trans_mut)
write.table(Trans,"Trans_gene_proi.txt",row.names=F,col.names=F,sep="\t",quote=F)
###cis gene priorization
Cis_mut=Per.M.AS$Cis
write.table(Cis_mut,"Cis_mut_AS.txt",row.names=F,col.names=F,sep="\t",quote=F)
Cis=Cis_mut[,2:3]
write.table(Cis,"Cis_gene_proi.txt",row.names=F,col.names=F,sep="\t",quote=F)

##################get the information of mutations and AS genes
trans<-read.csv("Trans_gene_proi.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
trans_AS<-read.csv("Trans_mut_AS.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
trans_sample<-merge(trans_AS,trans,by.x = "V2",by.y = "V1")
mutation<-read.csv("CHOL_sample_mut_code.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
colnames(trans_sample)<-c("mutgene","sample","ASgene")
colnames(mutation)<-c("sample","sampleID2","mutgene","ID","chr","start","end","type")
trans_mutation<-merge(mutation,trans_sample,by.x = c("mutgene","sample"),by.y = c("mutgene","sample"))
AS_out<-read.csv("Person_AS_event.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
trans_mutation_AS<-merge(trans_mutation,AS_out,by.x = c("sample","ASgene"),by.y = c("V1","V2"))
write.table(trans_mutation_AS,"DrAS_mut_AS_trans.txt",row.names=F,col.names=T,sep="\t",quote=F)

######cis pairs
cis<-read.csv("Cis_gene_proi.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
cis_AS<-read.csv("Cis_mut_AS.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
colnames(cis_AS)<-c("sample","mutgene","ASgene")
cis_sample<-merge(cis_AS,cis,by.x = "mutgene",by.y = "V1")
cis_mutation<-merge(mutation,cis_sample,by.x = c("mutgene","sample"),by.y = c("mutgene","sample"))
cis_mutation_AS<-merge(cis_mutation,AS_out,by.x = c("sample","ASgene"), by.y = c("V1","V2"))
write.table(cis_mutation_AS,"DrAS_mut_AS_cis.txt",row.names=F,col.names=T,sep="\t",quote=F)
