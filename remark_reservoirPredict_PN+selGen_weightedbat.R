"""
Babayan, Orton & Streicker
Predicting Reservoir Hosts and Arthropod Vectors from Evolutionary Signatures in RNA Virus Genomes
-- Reservoir host prediction from selected genomic features and phylogenetic neighborhoods

This script is used to generate training data to train nn models and also to generate nn input for prediction
"""

rm(list=ls())
setwd("") # Set local working directory where files are located

library(plyr)
library(h2o) # https://www.h2o.ai/products/h2o/
library(dplyr)
library(reshape2)
library(ape)
library(seqinr)
library(matrixStats)
library(hash) #used for weighted support score computation 
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

# Start h2o JVM
localh20<-h2o.init(nthreads = -1)  # Start a local H2O cluster using nthreads = num available cores

# define input output
queryVirusFeature = "SARS_related_CoV_representative_virusdata.csv"# "EPI_ISL_402131_VirusData.csv" # "Bat_Science_VirusData.csv" ##"Cov19VirusData.csv" # "Bat_VirusData.csv" to generate training data for nn model
queryVirusFasta = "SARS_related_CoV_representative_viruses.fasta" #"EPI_ISL_402131_sequences.fasta" # "Bat_Science_sequences.fasta" ## "Cov19Virus_seq.fasta" # "Bat_sequences.fasta" ##  "Bat_sequences.fasta" to generate training data for nn model
outputData = "SARS_related_CoV_representative_pseudoinput.csv" #"EPI_ISL_402131_pseudoinput.csv" # "dnn_batvirus_Science_pseudoinput.csv" ## "dnn_cov19virus_pseudoinput.csv" # "dnn_batvirus_pseudoinput.csv" ##"dnn_batvirus_input.csv" 
outputName = "SARS_related_CoV_representative_pseudonames.csv" #"EPI_ISL_402131_pseudonames.csv" # "dnn_batvirus_Science_pseudonames.csv" ##"dnn_cov19virus_pseudonames.csv" #"dnn_batvirus_pseudonames.csv" ## "dnn_batvirus_names.csv" 

# Read data from file
f1<-read.csv(file="BatPseudo_VirusData.csv",header=T) 
allP<-read.fasta(file ="BatPseudo_sequences.fasta",  seqtype = "DNA", as.string = TRUE, seqonly = F, strip.desc = T)
fis<-read.csv(file="featureImportance_reservoir.csv",header=T) 

# Feature definition
dinucs<-grep("[A|T|G|C|U]p[A|T|G|C|U]",names(f1),value=T)
cps<-grep(".[A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|X|Y]..[A|T|G|C|U]",names(f1),value=T)
aa.codon.bias<-grep(".Bias",names(f1),value=T)

# Feature selection (simplify dataset to required columns)
nfeats<-50
totalfeats<-length(fis$vimean)
f<-seq(from = totalfeats-(nfeats-1),to = totalfeats, by=1)
gen.feats<-as.character(fis$X[f]) # 50 most important target features as genome teatures gen.feats. note that decreasing=FALSE in order() in featureSelection.R
f1<-f1[,c("Virus.name","Genbank.accession","Reservoir","Viral.group","Vector.borne","Vector",gen.feats)] # extract target columns

# Remove orphans
f2<-subset(f1,f1$Reservoir!="Orphan")
f<-droplevels(f2)

# Group selection based on thresholds
t<-1 #15 # threshold for minimum sample size of groups 
s<-1 #.7 # proportion in the training set 
host.counts<-table(f$Reservoir)
min.t<-host.counts[host.counts>=t] # minimum number of viruses per host group
f_st3<-f[f$Reservoir %in% c(names(min.t)),]
f_st3<-droplevels(f_st3)
f_st3$SeqName2<-do.call(rbind,strsplit(as.character(f_st3$Genbank.accession),"[.]"))[,1]

# Number and names of host taxa
ntax<-length(unique(f_st3$Reservoir)) # number
bp<-as.character(sort(unique(f_st3$Reservoir))) # names
  
# Rare hosts. as target data whose blast result need to be computed. Bat_VirusData.csv
rare<-read.csv(file=queryVirusFeature,header=T) #rare<-read.csv(file="Cov19VirusData.csv",header=T) 
#rare<-f[!f$Reservoir %in% c(names(min.t)),]
rare<-droplevels(rare)
rare$SeqName2<-do.call(rbind,strsplit(as.character(rare$Genbank.accession),"[.]"))[,1]


# Train many models
set.seed(78910)
myhash <- hash(c("PTEROPODIDAE","VESPERTILIONIDAE"),c(1.0,1.0)) 
#myhash <- hash(c("PTEROPODIDAE","RarePterobats","RareVespbats","VESPERTILIONIDAE"),c(1.0,2.0,1.38,1.13)) 
#myhash <- hash(c("MOLOSSIDAE","PHYLLOSTOMIDAE","PTEROPODIDAE","RarePterobats","RareVespbats","VESPERTILIONIDAE"),c(4.0,3.27,1.0,2.0,5.14,1.13))
nloops<- 1 #550 
lr<-c()
md<-c()
sr<-c()
csr<-c()
nt<-c()
mr<-c()
accuracy.st3<-c()

pc.accuracy<-matrix(nrow=nloops,ncol=ntax)
#test.record<-matrix(nrow=ntest,ncol=nloops)
nfeatures<-length(gen.feats)+ntax
vimps<-matrix(nrow=nfeatures,ncol=nloops)

for (i in 1:nloops){
    # Stratified random. our model is meant to be trainDB independent
    trains<-f_st3 %>% group_by(Reservoir) %>% filter(Genbank.accession %in% sample(unique(Genbank.accession), ceiling(s*length(unique(Genbank.accession)))))
    trainSeqs<-allP[c(which(names(allP) %in% trains$Genbank.accession))] # pick sequences in the training set
    write.fasta(trainSeqs, names(trainSeqs), file.out="trainDB.fasta", open = "w", nbchar = 100, as.string = T)    
    # BLAST. should be repeated to make it not biased and in line with the bagging performance
    system("makeblastdb -in trainDB.fasta -dbtype nucl -parse_seqids -out allTrainingDB",intern=F)
    
    # Rare blast (take top 5 hits). queryVirusFasta. Bat_sequences.fasta. Bat_Science_sequences.fasta
    system("blastn -db allTrainingDB -query SARS_related_CoV_representative_viruses.fasta -out rareOut.out -num_threads 4 -outfmt 10 -max_target_seqs=5 -max_hsps 1 -reward 2 -task blastn -evalue 10 -word_size 8 -gapopen 2 -gapextend 2",inter=F,wait=FALSE)
    ### system("blastn -db allTrainingDB -query Cov19Virus_seq.fasta -out rareOut.out -num_threads 4 -outfmt 10 -max_target_seqs=5 -max_hsps 1 -reward 2 -task blastn -evalue 10 -word_size 8 -gapopen 2 -gapextend 2",inter=F,wait=FALSE)
    Sys.sleep(1) 
    # Summarize blast hits from rare virus set. This is what we need if we want to prepare input for rare set prediction.
    rBlast<-read.csv(file="rareOut.out",col.names = c("query acc.", "subject acc.", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end"," s. start"," s. end"," evalue"," bit score"),header=F)
    nvir<-length(unique(rBlast$query.acc.)) # read blast output rareout.out and get the number of query viruses in this data set
    virnames<-unique(rBlast$query.acc.) # names of query viruses
    ecutoff<-1E-3
    j=1
    d<-subset(rBlast,rBlast$query.acc.==virnames[j]) # pick one of the query viruses from blast result based on index j
    d2<-subset(d,d$X..identity<100) # exclude itself when doing pylogenetic analysis
    d2<-subset(d2,d2$X.evalue<ecutoff) # filter by e value
    # Assign equal probability across all hosts if there is no good blast hit
    for (z in 1:1){ # why do we need to seperate the 1st query virus from other query viruses? just for appending one by one?
      if (nrow(d2)==0){ # no hit at all. number of row in blast result is zero
        blast.uc<-rep(1/ntax,ntax) # number of host taxa
        blast.uc<-data.frame(t(blast.uc))
        colnames(blast.uc)<-sort(unique(trains$Reservoir)) # trains is replaced by f_st3 
        id<-as.character(virnames[j]) # j == 1
        blast.uc<-cbind(id,blast.uc)}
      else {
        dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F) # trains is replaced by f_st3 
        ##############################
        for (k in 1:length(dhost$X..identity)){
          dhost$pseudoidentity[k] <- dhost$X..identity[k] * myhash[[toString(dhost$Reservoir[k])]]
        }
        dhost$rel.support<-dhost$pseudoidentity/sum(dhost$pseudoidentity)
        ##############################
        #dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity)
        hosts<-tapply(dhost$rel.support,dhost$Reservoir,sum,na.rm=F) # result is like "NA NA 0.4111765 NA NA 0.5888235"
        hosts[is.na(hosts)]<-0
        hosts<-t(data.frame(hosts))
        hosts<-data.frame(hosts)
        id<-as.character(virnames[j]) # j == 1
        blast.uc<-cbind(id,hosts)}}
    
    if (nvir > 1){ # bug. it's been proven that if() is needed for the special case where nvir == 1 to avoid execution
    for (j in 2:nvir){ # if number of query viruses is 1, 2:1 will not execute anything?
      d<-subset(rBlast,rBlast$query.acc.==virnames[j]) # j starts from 2. pick one virus based on index j until all query viruses are done
      d2<-subset(d,d$X..identity<100) # exclude itself when doing pylogenetic analysis
      d2<-subset(d2,d2$X.evalue<ecutoff) # filter by e value
      if (nrow(d2)==0){ # the jth virus has no hit
        blast.uc.s<-rep(1/ntax,ntax) # assign an average probability to each host taxon
        blast.uc.s<-data.frame(t(blast.uc.s)) # convert to data frame
        colnames(blast.uc.s)<-sort(unique(trains$Reservoir)) # trains is replaced by f_st3 
        id<-as.character(virnames[j]) # query virus name
        blast.uc.s<-cbind(id,blast.uc.s) } # query virus name, probability of each host taxon based on blast result
      else {
        dhost<-merge(d2,trains,by.x="subject.acc.",by.y="Genbank.accession",all.x = T,all.y = F) # trains is replaced by f_st3 
        ##############################
        for (k in 1:length(dhost$X..identity)){
          dhost$pseudoidentity[k] <- dhost$X..identity[k] * myhash[[toString(dhost$Reservoir[k])]] # dhost$Reservoir[k] is a factor so need to convert to string
        }
        dhost$rel.support<-dhost$pseudoidentity/sum(dhost$pseudoidentity)
        ##############################
        #dhost$rel.support<-dhost$X..identity/sum(dhost$X..identity) # data for host dhost is the query's blast result that shows subject viruses, each subject virus' 50 traits
        hosts<-tapply(dhost$rel.support,dhost$Reservoir,sum,na.rm=T) # the reservori column as reference to operate in a group way, then sum up support score (percentage) for each host reservoir group (e.g on group basis), remove na value from result 
        hosts[is.na(hosts)]<-0 # replace na with 0
        hosts<-t(data.frame(hosts)) # convert to data frame then `t()` transposes the rows and columns of matrices
        hosts<-data.frame(hosts) # convert transposed result to data frame
        id<-as.character(d$query.acc.[1]) # id of query virus. things like NC_045512.2 NC_045512.2 NC_045512.2 NC_045512.2 NC_045512.2. so [1] is NC_045512.2
        blast.uc.s<-cbind(id,hosts)} # query virus, blast support score for each host
      blast.uc<-rbind(blast.uc,blast.uc.s)}} # bind by row, combine with result for the previous query virus 

    f1_rare<-merge(rare,blast.uc,by.x="Genbank.accession",by.y="id",all.x=T,all.y=T,sort=F) # query virus traits, blast result of query virus that shows support score for each host taxa. "accession" must match "id" when merging if all.x all.y both are "T"?
    set<-c(gen.feats,bp) # from merged info, only pick query virus' 50 features and all host taxa columns using taxa name 'bp'
    f1_rare<-f1_rare[,c(set)] # this is the input (50 features + scores for each host taxon) we need for rare set prediction


    # Convert to h2o data frames
    rar<-as.h2o(f1_rare)
    write.csv(f1_rare, file=outputData, row.names = F) #write.csv(f1_rare, file="dnn_cov19virus_input.csv", row.names = F) 
    write.csv(virnames, file=outputName, row.names = F) #write.csv(virnames, file="dnn_cov19virus_names.csv", row.names = F)

    # Clean up
    #rm(f1_rare)
    fn <- "rareOut.out"
    file.remove(fn) # have to remove otherwise by default blast does not rewrite the existing file

    # Identity the response column
    #y <- "Reservoir"

    # Identify the predictor columns
    #x <- setdiff(names(train), y)

    #load model according to model name and predict reservoir
    ##############################model_path <- "C:\\Users\\Bill\\Documents\\Rproject\\ViralHostPredictor\\model3rare\\bird_bat\\40"
    ##############################best_gbm <- h2o.loadModel(model_path)
    ##############################rar.pred <- h2o.predict(best_gbm, rar)
    
    # Retreive feature importance
    ##############################vi <- h2o.varimp(best_gbm)
    ##############################data2  <- vi[order(vi[,1],decreasing=FALSE),]
    #vimps[,i]<-data2[,4]

    # Rare predictions
    #rar.pred <- h2o.predict(best_gbm, rar)
    #df<-rar.pred[,c(2:(ntax+1))]
    #df2<-as.data.frame(df)
    #row.names(df2)<-rare$Virus.name
    #write.csv(df2,file=paste("ST5_RareVirus",i,".csv",sep="_"))

    # Clean up
    ##############################h2o.rm("gbm_grid")
    rm(oBlast,testBlast,allBlast,rBlast,trainSeqs,testSeqs,optSeqs,gbm_grid,best_gbm,train,test,opt,df2,optims)
}

#accs<-data.frame(accuracy.st3,pc.accuracy,lr,sr,md,csr)
#colnames(accs)[2:(ntax+1)]<-row.names(cm)
#row.names(vimps)<-data2$variable

# Write results summaries
#write.csv(vimps,file="Reservoir_PN+SelGen50_FI.csv",row.names = T)
#write.csv(accs,file="Reservoir_PN+SelGen50_out.csv",row.names=F)
