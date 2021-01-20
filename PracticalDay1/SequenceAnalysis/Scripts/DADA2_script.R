

#####
# Master Dades Omiques - Universitat de Vic.
# 2020-21, 26-27-28 Jan 2021
# OTU/ASV picking Analysis

# This R code aims to be an example of how a DADA2 workflow works
# DADA2 is an amplicon sequence variant analysis R packages
# https://www.nature.com/articles/ismej2017119
# Full Documentation here: https://benjjneb.github.io/dada2/index.html
#####

#### Using these version, other versions may also work, not checked.
library(dada2); packageVersion("dada2") ## Dada2 v.1.18.0
library(ShortRead); packageVersion("ShortRead") ## ShortRead 1.48.0
library(phyloseq); packageVersion("phyloseq") ## phyloseq 1.34.0
require(digest)


#### Set Working directory to the root directory where you have your data code
code_path<-("~/Documents/Work/Development/MasterDadesOmiques/PracticalDay1/SequenceAnalysis")
code_path
setwd(code_path)

#### Download RawData from

#### Download Training Data from
RawDataPath <- paste(code_path,"/RawData/",sep="")
TrainingPath<-paste0(code_path,"/Training/")


# Sort ensures forward/reverse reads are in same order
# Note that we are using a sequence subset for all samples of 10k sequences to make it feasible
fnFs <- sort(list.files(RawDataPath, pattern="_R1_10k.fastq.gz"))
fnRs <- sort(list.files(RawDataPath, pattern="_R2_10k.fastq.gz"))

# Infer sample names from file names
sample.names <- sapply(strsplit(fnFs, "_S"), `[`, 1)
#
# ##2020-12-23 18:23:26 MNJ commented out to read from metadatafile
# #metadata <-sapply(json_data, `[[`, "metadata")
# metadata<-read.table(paste(data_path,"Metadata/metadata.csv",sep=""))

#DONE-MNJ 2019-03-26 17:10:13 Log "how many files/samples we are processing"
##DBModule$insertLogs(sprintf("Processing %d samples",length(sample.names)),job_id,1)
fnFs <- file.path(RawDataPath, fnFs)
fnRs <- file.path(RawDataPath, fnRs)
plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])
filt_path <- file.path(code_path, "DADA2/filtered") # Place filtered files in filtered/ subdirectory

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

##filterAndTrim
## This takes less than a minute on a macbookpro
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,250),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE,verbose=T)

### How many sequences are passing the quality Filters
head(out)

plotQualityProfile(filtFs[1:2])
plotQualityProfile(filtRs[1:2])

### From filtered reads, we iteratively learn an error profile
### Should take less than 5 minutes on macboookpro
errF <- learnErrors(filtFs, multithread=TRUE,nbases=50000000,randomize=T,verbose=T)
errR <- learnErrors(filtRs, multithread=TRUE,nbases=50000000,randomize=T,verbose=T)

###
plotErrors(errF, nominalQ=TRUE)


### Using the error model that we have learnt from the data
### We can infer which reads are a result of true existence or a results from error from a true existing sequence.
###
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

##DBModule$insertLogs("Removing sequencing errors",job_id,1)
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]

##DBModule$insertLogs("Merging denoised forward and reverse reads.",job_id,1)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)


##DBModule$insertLogs("Removing chimeras",job_id,1)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
seqtab.nochim.pooled <- removeBimeraDenovo(seqtab, method="pooled", multithread=TRUE, verbose=TRUE)
### The following takes too long.
#seqtab.nochim.persample <- removeBimeraDenovo(seqtab, method="per-sample", multithread=TRUE, verbose=TRUE)


dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN),rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged",  "nonchim")
rownames(track) <- sample.names
head(track)
### Choose DB

### AssignTaxonomy
### Will need to chose taxonomy database from Analysis Configuration - Silva default
### We'll talk about the different databases in class.
list.files(paste(code_path,"/Training/",sep=""))
taxa.silva<-readRDS("taxa.silva.rds")
taxa.silva<- assignTaxonomy(seqtab.nochim, paste(code_path,"/Training/silva_nr_v132_train_set.fa.gz",sep=""), multithread=TRUE)

### Assign species if possible
### Since exact sequences are assigned we can try to assign these sequences to species level in some cases (specific genus, not generalizable)
### Assignment to species highly depends on the database of choice.
taxa.silva <- addSpecies(taxa.silva, paste0(code_path,"/Training/silva_species_assignment_v132.fa.gz",""),verbose=T)


###
DFForSeq<-colnames(seqtab.nochim)
system("rm DADA2/DADA_seq.fas")

for(seq in DFForSeq){
  write(paste(">",digest(seq),"\n",seq,sep=""),file="DADA2/DADA_seq.fas",append=T)
}


#### Taxa rownames are now the sequence of the ASV itself. This is not manageable downstream
#### We transform sequences to shorter strings using md5-algorithm so we can keep track
rownames(taxa.silva)<-sapply(rownames(taxa.silva),digest,algo="md5")
colnames(seqtab.nochim)<-sapply(colnames(seqtab.nochim),digest,algo="md5")


### We have the sequences for each ASV in DADA_seq.fas
### We can use it in any othter downstream analysis, but right now we want
### to create a phylogenetic tree to attach it to our data. This tree will allow
### calculation of different classes of phylogneetic-based distances (Unifrac, Faith's index,...)
### We build the alignment with mafft, mafft should be in your path (conda install -c bioconda mafft)
### Need to know where mafft executable is (typically in ~/miniconda3/bin or similar)
system("~/miniconda3/bin/mafft --auto --thread -1 DADA2/DADA_seq.fas > DADA2/DADA_seq.afa")

### Based on the alignmnent (afa file) we calculate the tree
### We use FastTree for that, using a generalized time-reversible model
### need to know where fastree executable is (conda install -c bioconda fasttree)
system("~/miniconda3/bin/fasttree -nt -gtr DADA2/DADA_seq.afa > DADA2/DADA_seq.nwk")

### Now we read the phylogenetic tree into a phylogenetic object (ape package)
treeSilva<-read_tree("DADA2/DADA_seq.nwk")

### Remove huge objects.
rm(derepFs)
rm(derepRs)


### How did we do for our mock positive control samples?
unqs.mock <- seqtab.nochim["Sample-Cpos1",]
unqs.mock <- seqtab.nochim["Sample-Cpos2",]
unqs.mock <- seqtab.nochim["Sample-Cpos3-1rep",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
rownames(unqs.mock)


#### Read metadata
metadata<-read.csv(paste0(code_path,"/Metadata/metadata.csv"))
rownames(metadata)<-metadata$SampleID

#### Create PhyloSeq Object
#### File named metadata.csv must exist and have SampleID variable for sample name
ps_silva<-phyloseq(otu_table(seqtab.nochim, taxa_are_rows=F),tax_table(taxa.silva),treeSilva)


### Clean phyloseq object from mock controls
sample_data(ps_silva)<-metadata
### Clean workspace
rm(DFForSeq,filtFs,filtRs,fnRs,fnFs,seq,track,mergers,seqtab.nochim.bkp,seqtab.nochim.pooled,taxa.silva.bkp)
# #asiggn metadata
# ##DBModule$insertLogs("Assign metadata if exists",job_id,1)
# if(!any(sapply(metadata,function(list) length(list)==0))) {
#   sampledata<-metModule$obtainSampleMet(metadata,sample.names)
#   otu_samples<-rownames(otu_table(ps_silva))
#   meta_samples<-rownames(sampledata)
#   if (all(length(otu_samples)==length(meta_samples)) && all(otu_samples==meta_samples)) {
#     sample_data(ps_silva)<-sampledata
#     sample_names(ps_silva)
#   }
# }
# ps_silva

### We keep the session
save.image(file=paste0(code_path,"/DADA2/DADA2_Rsession.RData"))

