##### Metagenomics course, 2018-19
##### Marc Noguera-Julian, PhD, IrsiCaixa, UVic AIDS Chair
##### 2018, Feb 19th.

##### In this part of the workshop we will analyze an otu table and the associated metadata
##### Data is obtained from a study realized in IrsiCaixa regarding the role of gut microbiome
##### in HIV infection and clinical evolution

##### We have metagenomic-16s data from 110 samples including healthy and HIV-infected patients

# Part 0. We need to setup the environment
#######   read the data and adapt metadata for further analysis

### Set working directory where this file is located
myWD<-NA
if(is.na(myWD)){
  print("Using rstudio hardcoded WD variable")
}
setwd(myWD)
#######   We will use the following packages: phyloseq, vegan, ggplot2 and DeSeq2
suppressPackageStartupMessages(library(phyloseq))### Bioconductors
suppressPackageStartupMessages(library(vegan))### Bioconductors
suppressPackageStartupMessages(library(ggplot2))### CRAN
suppressPackageStartupMessages(library(DESeq2)) ### Bioconductors
suppressPackageStartupMessages(require(doBy))### CRAN
suppressPackageStartupMessages(require(RColorBrewer))### CRAN
#######   gdata packages allows for microsofft office filetype read-in
suppressPackageStartupMessages(library(gdata))### CRAN
packageVersion("phyloseq")
packageVersion("vegan")
packageVersion("ggplot2")
packageVersion("DESeq2")
packageVersion("RColorBrewer")
packageVersion("doBy")


###################### Read OTU Data ##########################
### Since we are using 16s amplicon design it will be possible to
### taxonomically classify down to the genus level,
AvailableRanks<-c("Kingdom","Phylum","Class","Order","Family","Genus")

### Define the minimum number of OTU counts for a sample to be further analyzed
minSampleCountsB<-1000
### Read the metadata
metadataB<-read.csv("Metadata_Workshop_v1.0.csv")

### Let's inspect the metadata for starters
dim(metadataB)
summary(metadataB)
#minOtuCounts<-100

### We will do some processing of some variables to facilitate further analysis
rownames(metadataB)<-metadataB$SampleID
metadataB$RiskGroup<-metadataB$RiskGroup2
metadataB$gender<-metadataB$Gender
metadataB$HIVStatus<-as.character(metadataB$HIV_Status)
metadataB[metadataB$HIV_Status=="positive",]$HIVStatus<-"HIVpos"
metadataB[metadataB$HIV_Status=="negative",]$HIVStatus<-"HIVneg"
metadataB$HIVStatus<-as.factor(metadataB$HIVStatus)
metadataB$SexualPractice<-"HET"
metadataB[metadataB$RiskGroup2=="msm",]$SexualPractice<-"MSM"
metadataB$SexualPractice<-as.factor(metadataB$SexualPractice)
metadataB$City<-"Barcelona"

### Now we will import data derived from Mothur. For this, the minimal information is
### an OTU table an a taxonomical classification for each of the OTUS
### In addition we can also import a phylogenetic tree containing evolutionary information
### between all OTUs. We will use a cut-off of 97% similarity to define an OTU
### in agreement with parameters used for Mothur pipeline
xB<-import_mothur(mothur_constaxonomy_file="final.an.unique_list.0.03.cons.taxonomy",
                  mothur_shared_file="SampleRenamed.shared",
                  mothur_tree_file="OtuRepSeq.tree",
                  cutoff=0.03)
xB
### Attach the metadata
sample_data(xB)<-metadataB
xB
### Define available ranks in taxonomy file
colnames(tax_table(xB))<-AvailableRanks

### Filter-off samples that are below 1000 counts (minSampleCountsB)
barplot(colSums(otu_table(xB)),las=2,cex.names=0.6,main="#Counts/Sample")
xB<-prune_samples(sample_sums(xB)>minSampleCountsB,xB)
xB

### Let's plot the number of counts for the remaining samples
barplot(colSums(otu_table(xB)),las=2,cex.names=0.6,main="#Counts/Sample")

### We can se that we filtered only 1 sample
# otu_table()   OTU Table:         [ 16566 taxa and 109 samples ]
# sample_data() Sample Data:       [ 109 samples by 108 sample variables ]
# tax_table()   Taxonomy Table:    [ 16566 taxa by 6 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 16566 tips and 16564 internal nodes ]
### Note that we have a total of 16566 different OTUs
### This high level of diversity is partially caused by sequencing error

### End Read Barcelona Data #######################




### Coverage barplots to a pdf file, for the record
pdf("CoverageBarplots.pdf",paper="A4r")
barplot(colSums(otu_table(xB)),las=2,cex.names=0.6,main="#Counts/Sample")
dev.off()


#### Start Alpha - Diversity Analysis
### We will characterize alpha-diversity indices using a rarefied subset of 5000 counts
### in order to balance sampling between samples
### First we try to apply a prevalence filter
x.2.0 <- xB
wh0=genefilter_sample(x.2.0,filterfun_sample(function(x) x>1), A=0.01*nsamples(x.2.0))
x.2.0<-prune_taxa(wh0,x.2.0)

### We can now study some rarefaction analysis
### we will use vegan package functions but also phyloseq objects

### We can draw a rarefaction curve using otu_table from x.2.0 5000 subset object
sample_data(x.2.0)$HIVStatus
rarecurve(data.frame(t(otu_table(x.2.0))),step=20,col=sample_data(x.2.0)$HIVStatus,label=F,legend=T)
rarecurve(data.frame(t(otu_table(x.2.0))),step=20,col=sample_data(x.2.0)$gender,label=F,legend=T)
rarecurve(data.frame(t(otu_table(x.2.0))),step=20,col=sample_data(x.2.0)$SexualPractice,label=F,legend=T)

### Let's see what would happen if we had not applied a simple filter
x.2.0 <- xB
rarecurve(data.frame(t(otu_table(x.2.0))),step=20,col=sample_data(x.2.0)$HIVStatus,label=F,legend=T)

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(x.2.0),
               MARGIN = ifelse(taxa_are_rows(x.2.0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(x.2.0),
                    tax_table(x.2.0))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(x.2.0, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(x.2.0),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
ggsave("Abundance_vs_Prevalence_Phylum.pdf")

### Let's go back
x.2.0 <- xB

### We randomly subset 5000 counts from each sample and discard all samples with fewer counts
x.2.0<-rarefy_even_depth(x.2.0,3000)  ##### RRarefy to 5000 counts

##########Generate Diversity Plots
p<-plot_richness(x.2.0,"HIVStatus",measures=c("Shannon","Simpson","InvSimpson"))
p+geom_boxplot(aes(fill=HIVStatus))
ggsave("DiversityByHIVstatus.pdf")
p<-plot_richness(x.2.0,"RiskGroup",measures=c("Shannon","Simpson","InvSimpson"))
p+geom_boxplot(aes(fill=RiskGroup))
ggsave("DiversityByRiskGroup.pdf")
p<-plot_richness(x.2.0,"SexualPractice",measures=c("Shannon","Simpson","InvSimpson"))
p+geom_boxplot(aes(fill=SexualPractice))
ggsave("DiversityBySexualPractice.pdf")
p<-plot_richness(x.2.0,"gender",measures=c("Shannon","Simpson","InvSimpson"))
p+geom_boxplot(aes(fill=gender))
ggsave("DiversityBygender.pdf")

########## Generate Richness Plots
p<-plot_richness(x.2.0,"HIVStatus",measures=c("Observed","Chao1","ACE"))
p+geom_boxplot(aes(fill=HIVStatus))
ggsave("RichnessByHIVstatus.pdf")
p<-plot_richness(x.2.0,"RiskGroup",measures=c("Observed","Chao1","ACE"))
p+geom_boxplot(aes(fill=RiskGroup))
ggsave("RichnessByRiskGroup.pdf")
p<-plot_richness(x.2.0,"SexualPractice",measures=c("Observed","Chao1","ACE"))
p+geom_boxplot(aes(fill=SexualPractice))
ggsave("RichnessBySexualPractice.pdf")
p<-plot_richness(x.2.0,"gender",measures=c("Observed","Chao1","ACE"))
p+geom_boxplot(aes(fill=gender))
ggsave("RichnessBygender.pdf")

my.sampleData<-data.frame(sample_data(x.2.0))
my.sampleData
er.x.2.0<-estimate_richness(x.2.0)
rownames(er.x.2.0)<-rownames(my.sampleData)

er.x.2.0<-merge(my.sampleData,er.x.2.0,by="row.names",all.x=T)
rownames(er.x.2.0)<-rownames(my.sampleData)
sample_data(x.2.0)<-er.x.2.0

### Plots are nice and look good (underline some trends, but we need statistics)
### Generate Diversity Stats and capture them in a txt file
capture.output(file="DivAndRichnessStats.txt",paste())
for(covar in c("HIVStatus","RiskGroup","SexualPractice","gender")){
  for (measure in c("Observed","Chao1","ACE","Shannon","Simpson","InvSimpson")){
    capture.output(file="DivAndRichnessStats.txt",paste("############",measure,"by",covar),append=T)
    if(length(unique(er.x.2.0[,covar]))==2){
      myttest<-t.test(er.x.2.0[,measure]~er.x.2.0[,covar])
      capture.output(file="DivAndRichnessStats.txt",myttest,append=T)
    }else{
      myanova<-aov(er.x.2.0[,measure]~er.x.2.0[,covar])
      capture.output(file="DivAndRichnessStats.txt",myanova,append=T)
      capture.output(file="DivAndRichnessStats.txt",summary(myanova),append=T)
      capture.output(file="DivAndRichnessStats.txt",TukeyHSD(myanova),append=T)
    }
  }
}

########################################## End of Alpha-Diversity Analysis ##########################################
#####################################################################################################################

#####################################################################################################################
######################################### Taxonomical analysis ######################################################
### Let's make an initial taxonomy description of the data
### Cumulative Stacked Barplots at Phylum, Genus and Species
### 3.0
x.3.0<-xB
### We apply a more stringent filter, basically we are not interested in rare OTUs but on
### general main trends of taxonomical composition
wh0=genefilter_sample(x.3.0,filterfun_sample(function(x) x>5), A=0.1*nsamples(x.3.0))
x.3.0<-prune_taxa(wh0,x.3.0)
### tax_glom function agglomerates/collapses all OTU belonging to the same taxonomical level
x.3.0.phylum<-tax_glom(x.3.0,taxrank="Phylum")
x.3.0.genus<-tax_glom(x.3.0,taxrank="Genus")
x.3.0.genus<-subset_taxa(x.3.0.genus,Genus != "unclassified")
x.3.0.genus<-subset_taxa(x.3.0.genus,Genus != "Incertae_Sedis")
### Convert to data.frame for easier manipulation
ps.melt.x.3.0.phylum<-psmelt(x.3.0.phylum)
ps.melt.x.3.0.genus<-psmelt(x.3.0.genus)

### Let's look at the phylum data.frame
summary(ps.melt.x.3.0.phylum)

### We will calculate relative abundances as proportions for further analysis
ps.melt.x.3.0.phylum$AbundanceProportion <- ave(ps.melt.x.3.0.phylum$Abundance,list(ps.melt.x.3.0.phylum[,"SampleID"]), FUN=function(L) L/sum(L))
ps.melt.x.3.0.genus$AbundanceProportion <- ave(ps.melt.x.3.0.genus$Abundance,list(ps.melt.x.3.0.genus[,"SampleID"]), FUN=function(L) L/sum(L))

reverse.levels <- function(x) {
  if(is.factor(x)) {
    x <- factor(as.character(x), levels=rev(levels(x)), ordered=TRUE)
  } else if(is.data.frame(x)) {
    for(i in seq_along(x)) {
      if(is.factor(x[,i])) {
        x[,i] <- factor(as.character(x[,i]), levels=rev(levels(x[,i])), ordered=TRUE)
      } else {
        warning(paste0('Column ', i, ' is not a factor.'))
      }
    }
  } else {
    stop(paste0('Unsupported format: ', class(x)))
  }
  return(x)
}

### Barplots at phylum. We use decreasing phylum abundances as X order
### and the phylum abundance of first sample as y order
### This careful use of X/Y order will already reveal patterns
### Define Level Order for X axis (SampleID)
my.levels=orderBy(~-AbundanceProportion,data=ps.melt.x.3.0.phylum[ps.melt.x.3.0.phylum$Phylum=="Bacteroidetes",])$SampleID
ps.melt.x.3.0.phylum$SampleID<-factor(ps.melt.x.3.0.phylum$SampleID,
                                      levels=my.levels,ordered=T)

### Define Level Order for Y axis (phylum,genus,species)
my.levels<-orderBy(~-AbundanceProportion,data=ps.melt.x.3.0.phylum[ps.melt.x.3.0.phylum$SampleID == as.character(my.levels[1]),])$Phylum
ps.melt.x.3.0.phylum$Phylum<-factor(ps.melt.x.3.0.phylum$Phylum,levels=my.levels,ordered=T)

### We also take care of colors
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colorCount=length(unique(ps.melt.x.3.0.phylum$Phylum))
getPalette=colorRampPalette(brewer.pal(12,"Set3"))
p<-ggplot(ps.melt.x.3.0.phylum,aes(x=SampleID,y=AbundanceProportion,order=ps.melt.x.3.0.phylum$Phylum))
p+geom_bar(stat="identity",aes(fill=Phylum))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+labs(title="Phylum Level")+
  scale_fill_manual(values=getPalette(colorCount))
### Let's save the plot
ggsave("x.3.0.phylum.barplot.ordered.pdf")

#Barplots at genus
# Define Level Order for X axis (SampleID)
my.levels=orderBy(~-AbundanceProportion,data=ps.melt.x.3.0.genus[ps.melt.x.3.0.genus$Genus=="Bacteroides",])$SampleID
ps.melt.x.3.0.genus$SampleID<-factor(ps.melt.x.3.0.genus$SampleID,
                                      levels=my.levels,ordered=T)
#Define Level Order for Y axis (genus,genus,species)
my.levels<-orderBy(~-AbundanceProportion,data=ps.melt.x.3.0.genus[ps.melt.x.3.0.genus$SampleID == as.character(my.levels[1]),])$Genus
ps.melt.x.3.0.genus$Genus<-factor(ps.melt.x.3.0.genus$Genus,levels=my.levels,ordered=T)

ps.melt.x.3.0.genus<-subset(ps.melt.x.3.0.genus,ps.melt.x.3.0.genus$AbundanceProportion>0.02)
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colorCount=length(unique(ps.melt.x.3.0.genus$Genus))
getPalette=colorRampPalette(brewer.pal(12,"Set3"))
p<-ggplot(ps.melt.x.3.0.genus,aes(x=SampleID,y=AbundanceProportion,order=Genus))
p+geom_bar(stat="identity",aes(fill=Genus))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+labs(title="genus Level")+
  scale_fill_manual(values=getPalette(colorCount))+guides(fill=guide_legend(ncol=2))
ggsave("x.3.0.genus.barplot.ordered.pdf")

### We see that there is a trend of decreasing abundance of certain genus across all the sample set
### however we do not really know now what this trend is related to


### We can assess Differential Abundance, related to any variable in our metadata
### 3.0
x.3.0<-xB
wh0=genefilter_sample(x.3.0,filterfun_sample(function(x) x>1), A=0.1*nsamples(x.3.0))
x.3.0<-prune_taxa(wh0,x.3.0)
### tax_glom function agglomerates/collapses all OTU belonging to the same taxonomical level
x.3.0.phylum<-tax_glom(x.3.0,taxrank="Phylum")
x.3.0.genus<-tax_glom(x.3.0,taxrank="Genus")

ps.melt.x.3.0.phylum<-psmelt(x.3.0.phylum)
ps.melt.x.3.0.genus<-psmelt(x.3.0.genus)

### Let's transform our data in different ways
ps.melt.x.3.0.phylum$Abundance<-floor(ps.melt.x.3.0.phylum$Abundance)
### Add a pseudocount to allow for log transformation
ps.melt.x.3.0.phylum$Abundance<-ps.melt.x.3.0.phylum$Abundance+1
ps.melt.x.3.0.phylum$AbundanceProportion <- ave(ps.melt.x.3.0.phylum$Abundance,list(ps.melt.x.3.0.phylum[,"SampleID"]), FUN=function(L) L/sum(L))
ps.melt.x.3.0.phylum$AbundanceProportionLog <- log10(ps.melt.x.3.0.phylum$AbundanceProportion)
ps.melt.x.3.0.phylum$AbundanceNorm<-floor(10000*ps.melt.x.3.0.phylum$AbundanceProportion)

ggplot(ps.melt.x.3.0.phylum,aes(x=Phylum,y=AbundanceProportionLog,fill=Phylum))+geom_boxplot()+facet_wrap(~HIVStatus)+ theme(axis.text.x = element_text(angle = 90))
ggplot(ps.melt.x.3.0.phylum,aes(x=Phylum,y=AbundanceProportionLog,fill=Phylum))+geom_boxplot()+facet_wrap(~RiskGroup)+ theme(axis.text.x = element_text(angle = 90))
ggplot(ps.melt.x.3.0.phylum,aes(x=Phylum,y=AbundanceProportionLog,fill=HIVStatus))+geom_boxplot()+facet_wrap(~SexualPractice)+ theme(axis.text.x = element_text(angle = 90))
ggplot(ps.melt.x.3.0.phylum,aes(x=Phylum,y=AbundanceProportionLog,fill=Gender))+geom_boxplot()+facet_wrap(~HIVStatus)+ theme(axis.text.x = element_text(angle = 90))

ps.melt.x.3.0.genus$Abundance<-floor(ps.melt.x.3.0.genus$Abundance)
ps.melt.x.3.0.genus$Abundance<-ps.melt.x.3.0.genus$Abundance+1
ps.melt.x.3.0.genus$AbundanceProportion <- ave(ps.melt.x.3.0.genus$Abundance,list(ps.melt.x.3.0.genus[,"SampleID"]), FUN=function(L) L/sum(L))
ps.melt.x.3.0.genus$AbundanceProportionLog <- log10(ps.melt.x.3.0.genus$AbundanceProportion)
ps.melt.x.3.0.genus$AbundanceNorm<-floor(10000*ps.melt.x.3.0.genus$AbundanceProportion)

### To plot at Genus level we first need to select the most abundant genus
### plotting all of them is a terrible mess
levels(ps.melt.x.3.0.genus$Genus)
meanByGenus<-aggregate(ps.melt.x.3.0.genus$AbundanceProportion,list(ps.melt.x.3.0.genus$Genus),mean)
colnames(meanByGenus)<-c("Genus","meanAbundance")
meanByGenus<-meanByGenus[with(meanByGenus, order(-meanAbundance, Genus)), ]


#### Plot relative logProportions of most 10 abundant genus
ggplot(ps.melt.x.3.0.genus[ps.melt.x.3.0.genus$Genus %in% meanByGenus[1:10,1],],aes(x=Genus,y=AbundanceProportionLog,fill=Genus))+geom_boxplot()+facet_wrap(~HIVStatus)+ theme(axis.text.x = element_text(angle = 90))
ggplot(ps.melt.x.3.0.genus[ps.melt.x.3.0.genus$Genus %in% meanByGenus[1:10,1],],aes(x=Genus,y=AbundanceProportionLog,fill=Genus))+geom_boxplot()+facet_wrap(~RiskGroup)+ theme(axis.text.x = element_text(angle = 90))

########################################################################
############ Start Genus comparison barplots #########################
########################################################################
### We have seen that there are differences but we need to see which these differences are
###Barplots for Significant Genus difference
x.4.0<-xB
wh0=genefilter_sample(x.4.0,filterfun_sample(function(x) x>1), A=0.01*nsamples(x.4.0))
x.4.0<-prune_taxa(wh0,x.4.0)
x.4.0 = transform_sample_counts(x.4.0, function(x) ((x/sum(x))))
x.4.0<-tax_glom(x.4.0,taxrank="Genus")
psmelt.x.4.0.genus<-psmelt(x.4.0)

##### Barplots for significant genus related to Sexual Practice
mysignificantGenus<-vector()
my.p.values.vector<-vector()
# Find significant genus (p<0.01)
mainDir<-"./"
subDir<-"StatisticalTests"
if (file.exists(subDir)){
  #setwd(file.path(mainDir, subDir))
} else {
  dir.create(file.path(mainDir, subDir))
  #setwd(file.path(mainDir, subDir))
}
### For every Genus we will calculate a Wilcoxon Rank-based test between dicotomic variable SexualPractice
for(genus in as.vector(unique(psmelt.x.4.0.genus$Genus))){
  mytest<-wilcox.test(Abundance~SexualPractice,data=psmelt.x.4.0.genus[psmelt.x.4.0.genus$Genus==genus,])
  capture.output(mytest,file=paste(mainDir,"/",subDir,"/",genus,"_SexualPractice_Wilcoxon.txt",sep=""))
  my.p.values.vector<-c(my.p.values.vector,mytest$p.value)
}

### We need to correct for multiple test error, we use benjamin-hochberg
my.p.values.vector.adj<-p.adjust(my.p.values.vector,method="BH")
mysignificantGenus<-unique(psmelt.x.4.0.genus$Genus)[my.p.values.vector.adj<0.05]
ps.x.4.0.genus<-psmelt.x.4.0.genus[psmelt.x.4.0.genus$Genus %in% mysignificantGenus,]

p<-ggplot(ps.x.4.0.genus,aes(x=Genus,y=Abundance))

#Plot B&W boxplot
p+geom_boxplot(aes(fill=SexualPractice,stat="identity"),notch=F,position="dodge")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ scale_fill_grey(start = 0.6, end = 1)+
  theme_bw() + theme(legend.title=element_blank(),legend.text=element_text(size=22),legend.position=c(.5, .5))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size=0.9),
        axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(angle = 45, hjust = 1,size=8),axis.title=element_text(size=8,face="bold"))
ggsave("GenusBoxplot_bySexualPractice_bw.pdf")  ################### <----  Save Pdf
#Plot Color barplot
p+geom_boxplot(aes(fill=SexualPractice,stat="identity"),notch=F,position="dodge")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw() + theme(legend.title=element_blank(),legend.text=element_text(size=22),legend.position=c(.5, .5))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size=0.9),
        axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(angle = 45, hjust = 1,size=8),axis.title=element_text(size=8,face="bold"))+scale_y_continuous(trans='sqrt')
ggsave("GenusBoxplot_bySexualPractice_color.pdf")  ################### <----  Save Pdf

require("ggpubr")
p+geom_boxplot(aes(fill=SexualPractice),notch=F,position="dodge")+facet_wrap(~Genus,scales="free")+
  theme_bw()+
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16))+
  theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
  theme(legend.text=element_text(size=16,face="italic"),legend.title=element_text(size=16))+
  theme(plot.title=element_text(lineheight=1,face="bold",size=19))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size=0.9),
        axis.text.x = element_blank(),
        axis.text.y = element_text( hjust = 1,size=8),
        axis.title=element_text(size=12,face="bold"),
        strip.text = element_text(size=8,face="bold"))+
  scale_colour_manual(values=c("darkred","darkolivegreen","dodgerblue2"))+scale_y_continuous(trans='sqrt')
ggsave("GenusBoxplot_bySexualPractice_color_WrapByGenus.pdf")   ################### <----  Save Pdf


###
statisticValue<-function(value=NULL){
  if(value>0.1){
    return (paste("=",round(value,digits=2)))
  }
  else if(value>=0.05){
    return (paste("=",round(value,digits=2)))
  }
  else if(value>=0.01){
    return (paste("<0.05"))
  }
  else if(value>=0.001){
    return (paste("<0.01"))
  }
  else if(value>=0.0001){
    return(paste("<0.001"))
  }
  else{return (paste("<0.0001"))}
}
boxplotNumericByGroup<-function(mydata,category,variable,nbvariable,test,Rank=NULL){
  # if(is.null(Rank)){
  #    title<-paste(as.character(variable), " by ",category)
  #    fileForOutput<-paste(as.character(variable),"by",category,sep="_")
  #  }else{
  Phylum<-unique(as.vector(mydata[,Rank]))
  title<-paste(as.character(variable), " of\n ",Phylum,as.character(Rank)," by ",category)
  fileForOutput<-paste(as.character(variable),"of",Phylum,as.character(Rank),"by",category,sep="_")
  #  }
  numberOfLevels<-length(unique(mydata[,category]))
  colorsPlotRight<-c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999")
  require(gridExtra)
  if(numberOfLevels==2){
    test<-t.test(mydata[,variable]~mydata[,category])
    testString<-paste("Students t-test. p-value",statisticValue(test$p.value))
  }
  if(numberOfLevels>=3){
    require(MASS)
    glm.nb.model<-glm.nb(mydata[,nbvariable]~mydata[,category],method="glm.fit")
    glm.nb.aov<-aov(glm.nb.model)
    test<-aov(mydata[,variable]~mydata[,category])
    tukey.test<-TukeyHSD(test)
    print(tukey.test$mydata)
    mymatrix<-tukey.test$mydata
    mymatrix<-mymatrix[,c("diff","p adj")]
    mymatrix<-as.matrix(mymatrix)
    #mymatrix<-round(mymatrix,3)
    for(i in 1:nrow(mymatrix)){
      print(mymatrix[i,"p adj"])
      mymatrix[i,"p adj"]<-statisticValue(as.numeric(mymatrix[i,"p adj"]))
    }
    text2.df<-as.table(mymatrix)
    testString<-paste("ANOVA PR(>F)", statisticValue(summary(test)[[1]][["Pr(>F)"]][[1]]))
    testString<-paste(testString,"\n","NegBin ANOVA PR(>F)",statisticValue(summary(glm.nb.aov)[[1]][["Pr(>F)"]][[1]]))
  }
  #  if(test=="lm"){
  #    model<-lm(mydata[,variable]~mydata[,category])
  #  }
  mydata$xaxis<-mydata[,category]
  mydata$yvalue<-mydata[,variable]
  p<-ggplot(mydata,aes(x=xaxis,y=yvalue,fill=as.factor(xaxis)))+
    geom_boxplot()+geom_jitter(color="DarkRed")+
    #ggtitle(title)+
    xlab(category)+ylab(variable)+
    ylim(min(mydata$yvalue-1),1)+
    scale_fill_manual(values=colorsPlotRight[1:numberOfLevels])+
    theme(legend.position=c(1,1),legend.justification=c(1,1))+
    #annotate("text",x=numberOfLevels/2.5,y=0.5,label=testString,size=3)+
    annotate("text",x=1,y=0,label=testString,size=3)+
    annotate("text",y=0.4,x=1.5,label=title,size=4)
  plotRight<-ggplot(mydata,aes(yvalue,fill=xaxis))+geom_density(alpha=.5)+
    coord_flip()+scale_fill_manual(values=colorsPlotRight[1:numberOfLevels])+
    theme(legend.position="none")+
    xlim(min(mydata$yvalue-1),max(mydata$yvalue+1))
  if(numberOfLevels>=3){
    p<-p+annotation_custom(tableGrob(text2.df), ymin=min(mydata$yvalue)-1, ymax=min(mydata$yvalue), xmax=numberOfLevels/1.2, xmin=numberOfLevels/2)
  }  #p2<-tableGrob(text2.df)
  #grid.arrange(p2,p,main="prova",ncol=2)
  fileForPlot <- paste(fileForOutput,".pdf")
  pdf(fileForPlot,paper="a4r")
  grid.arrange(p,plotRight,nrow=1,ncol=2,widths=c(4,1),heights=c(4))
  #p2<-ggplot(p2)
  #ggsave(filename = fileForPlot,dpi=600,width=11, height=8.5)
  dev.off()
  grid.arrange(p,plotRight,nrow=1,ncol=2,widths=c(4,1),heights=c(4))
  #print(p2)
  #return(p2)
}
compareTwoGroups<-function(mydata=NULL,variable=NULL,category1=NULL,category2=NULL,fileForPlot=NULL,minCounts=500,maxAlpha=0.01, design=NULL){
  fileForPlot=paste("NegBin_DiffTest_",as.character(deparse(substitute(mydata))),"_",variable,"_",category1,"vs",category2,".pdf",sep="")
  print(paste("Output file for Plots: ",fileForPlot))
  fileForTable=paste("NegBin_DiffTest_",as.character(deparse(substitute(mydata))),"_",eval(variable),"_",eval(category1),"vs",eval(category2),".txt",sep="")
  print(paste("Output file for Table: ",fileForTable))
  stringForTitle=paste(variable,"/",category1,"vs",category2)
  kostic <- mydata
  LastRank<-rank_names(mydata)[length(rank_names(mydata))][[1]]
  require(phyloseq)
  require(DESeq2)
  kostic <- prune_samples(sample_sums(kostic) > minCounts, kostic)
  diagdds = phyloseq_to_deseq2(kostic, as.formula(design))
  #print(diagdds)
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(diagdds), 1, gm_mean)
  diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
  #colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("untreated","treated"))
  #colData(diagdds)$condition<-factor(colData(diagdds)$condition,levels=c(category1,category2))
  diagdds = DESeq(diagdds, fitType="parametric",test="Wald")
  res=results(diagdds,contrast=c(variable,category1,category2))
  #res=results(diagdds)
  print(res)
  res = res[order(res$padj, na.last=NA), ]
  sigtab = res[(res$padj < maxAlpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(kostic)[rownames(sigtab), ], "matrix"))
  sigtab$OtuID<-rownames(sigtab)
  head(sigtab)
  write.table(sigtab,file=fileForTable,sep="\t")
  #Cleanup for Positive enrichment in csigtabarcinoma
  posigtab=sigtab
  #posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
  #posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
  posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", rank_names(mydata))]
  library("ggplot2")
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  sigtabgen=sigtab

  #sigtabgen = subset(sigtab, !is.na(Genus))
  #sigtabgen = subset(sigtabgen, sigtabgen$Genus != "unclassified")
  # Phylum order
  x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
  # Genus order
  if(as.character(LastRank)== "Genus"){
    sigtabgen$LastRank<-sigtabgen[,"Genus"]
  }else{
    sigtabgen$LastRank<-paste(sigtabgen[,"Genus"]," ",sigtabgen[,as.character(LastRank)])
  }
  x = tapply(sigtabgen$log2FoldChange, sigtabgen$LastRank, function(x) max(x))
  x = sort(x, TRUE)

  sigtabgen$LastRank = factor(as.character(sigtabgen$LastRank), levels=names(x))
  sigtabgen$log2Counts<-log2(sigtabgen$baseMean)
  sigtabgen$alpha<- 1 - sigtabgen$padj
  #pdf(fileForPlot)
  p<-ggplot(sigtabgen,aes(x=LastRank,y=log2FoldChange))
  p<-p+geom_point(aes(colour=Phylum,size=log2Counts,alpha=alpha))
  p<-p+scale_size_continuous(range=c(1,20))
  #+geom_point(aes(size=sigtabgen$log2Counts))+scale_size_area()
  p<-p+theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=0.5,size=10))
  p<-p+theme(legend.key.size=unit(1,"cm"))
  p<-p+ ggtitle(paste(stringForTitle," Data:",as.character(deparse(substitute(mydata))))) +
    theme(plot.title = element_text(lineheight=.7, face="bold"))
  print(p)
  #ggplot(sigtabgen, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6)+scale_size(range=c(1,5))+
  # theme(axis.text.x = element_text(angle = -90, hjust = 0,size=3, vjust=0.5), legend.key.size=unit(0.5,"cm"),legend.text=element_text(size=3))
  ggsave(filename = fileForPlot,dpi=600,width=11, height=8.5)
  return(sigtab)
  #dev.off()
}


my.subset<-subset(ps.melt.x.3.0.phylum,ps.melt.x.3.0.phylum$Phylum=="Firmicutes")
p1<-boxplotNumericByGroup(my.subset,category="RiskGroup",variable="AbundanceProportionLog",nbvariable="Abundance",Rank="Phylum")
my.subset<-subset(ps.melt.x.3.0.genus,ps.melt.x.3.0.genus$Genus=="Prevotella")
p1<-boxplotNumericByGroup(my.subset,category="RiskGroup",variable="AbundanceProportionLog",nbvariable="Abundance",Rank="Genus")
p1<-boxplotNumericByGroup(my.subset,category="HIVStatus",variable="AbundanceProportionLog",nbvariable="Abundance",Rank="Genus")


### Let's move one step further and screen all genus for significant difference by using
### a more complex (and adequate) statistical framework using DESeq2 package and negative binomial
### distribution fits to detect over(under)-represented genus in a dichotomic condition

compareTwoGroups(mydata=x.3.0.genus,variable="HIVStatus",category1="HIVpos",category2="HIVneg",design=~HIVStatus,maxAlpha=0.01)
compareTwoGroups(mydata=x.3.0.genus,variable="SexualPractice",category1="MSM",category2="HET",design=~SexualPractice,maxAlpha=0.01)

#### End of Bacterial taxonomical analysis

### Start Ordination Analysis
### For this we need to use ecological distance. We will be using both Bray-Curtis distance
### and WUnifrac phylogenetic distance for this tutorial

#Let put our data aside
x.4.0<-xB

### Remember that we have more than 16k OTUs in our data. Most of them are unique to a single sample
### thus indicating that they are probably artifact. Usually, these OTU don't have a great impact on distance
### calculations for this particular reason

#####Intended for de novo OTU-picking only, we keep OTUs that appear at least in two different samples
wh0=genefilter_sample(x.4.0,filterfun_sample(function(x) x>0), A=0.01*nsamples(x.4.0))
x.4.0<-prune_taxa(wh0,x.4.0)
#x.4.0<-tax_glom(x.4.0,taxrank="Genus")

### We transform the data to proportion abundances
x.4.0 = transform_sample_counts(x.4.0, function(x) ((x/sum(x))))

### We will be using NMDS which used an initial random seed for iteration
### Since we want to be able to repeat the same analysis we fix a seed for analysis
set.seed(12345)


### NMDS Ordination with Bray-Curtis distance
x.4.0.ord<-ordinate(x.4.0,"NMDS",distance="bray",trymax=200)
capture.output(file="NMDS_Bray_proportions_ordinfo.txt",x.4.0.ord)
pdf("NMDS_Bray_proportions_stressplot.pdf")
stressplot(x.4.0.ord)
dev.off()
### Let's take a look at stressplot. Is it good?
stressplot(x.4.0.ord)
### Simple Plot
p.4.0.samples=plot_ordination(x.4.0,x.4.0.ord)
p.4.0.samples

### Let's map some metadata with coloured and ellipses
p.4.0.samples=plot_ordination(x.4.0,x.4.0.ord)
p.4.0.samples  + geom_point(aes(color=SexualPractice,shape=HIVStatus,size=Age)) +ggtitle("Unconstrained NMDS(Bray-Curtis)")+
  scale_colour_manual(values=c("darkred","darkolivegreen","dodgerblue2"))+
  theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),axis.text.x=element_text(size=16),axis.text.y=element_text(size=16))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  theme(legend.text=element_text(size=12,face="italic"),legend.title=element_text(size=12))+
  stat_ellipse(geom="polygon",alpha=0.25,aes(fill=SexualPractice),level=0.95)+
  theme(plot.title=element_text(lineheight=1,face="bold",size=15))
ggsave("NMDS_Bray_Proportion_Colour_withClusters_SizeByAge_BCN.pdf",) ################### <----  Save Pdf

### Now using HIV STatus
p.4.0.samples  + geom_point(size = 5.5,aes(color=HIVStatus,fill=HIVStatus,shape=SexualPractice)) +ggtitle("Unconstrained NMDS(Bray-Curtis)")+
  stat_ellipse(geom="polygon",alpha=0.25,aes(fill=HIVStatus),level=0.95)+
  scale_colour_manual(values=c("darkred","darkolivegreen","dodgerblue2"))+
  theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),axis.text.x=element_text(size=16),axis.text.y=element_text(size=16))+
  theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
  theme(legend.text=element_text(size=16,face="italic"),legend.title=element_text(size=16))+
  theme(plot.title=element_text(lineheight=1,face="bold",size=19))
ggsave("NMDS_Bray_Proportion_Colour_withClusters_byHIVStatus.pdf")

### Let's do a similar analysis but using phylogenetic OTU-OTU relationship
### using unifrac distances
x.4.0.ord<-ordinate(x.4.0,"NMDS",distance="wunifrac")
capture.output(file="NMDS_WUnifrac_proportions_ordinfo.txt",x.4.0.ord)
pdf("NMDS_WUnifrac_proportions_stressplot.pdf")
stressplot(x.4.0.ord)
dev.off()

### Simple Plot
p.4.0.samples=plot_ordination(x.4.0,x.4.0.ord)
p.4.0.samples

### Now with metadata mapping
p.4.0.samples=plot_ordination(x.4.0,x.4.0.ord)
p.4.0.samples  + geom_point(size = 5.5,aes(color=SexualPractice,fill=RiskGroup,shape=HIVStatus,size=Age)) +ggtitle("Unconstrained NMDS(Weighted Unifrac)")+
  scale_colour_manual(values=c("darkred","darkolivegreen","dodgerblue2"))+
  theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),axis.text.x=element_text(size=16),axis.text.y=element_text(size=16))+
  theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
  theme(legend.text=element_text(size=16,face="italic"),legend.title=element_text(size=16))+
  stat_ellipse(geom="polygon",alpha=0.25,aes(fill=SexualPractice),level=0.95)+
  ggtitle("Genus level")+
  theme(plot.title=element_text(lineheight=1,face="bold",size=19))
ggsave("NMDS_Wunifrac_Proportion_Colour_withClusters.pdf")

p.4.0.samples  + geom_point(size = 5.5,aes(color=HIVStatus,fill=HIVStatus,shape=SexualPractice,size=Age)) +ggtitle("Unconstrained NMDS(Weighted UniFrac)")+
  stat_ellipse(geom="polygon",alpha=0.25,aes(fill=HIVStatus),level=0.95)+
  scale_colour_manual(values=c("darkred","darkolivegreen","dodgerblue2"))+
  theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),axis.text.x=element_text(size=16),axis.text.y=element_text(size=16))+
  theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
  theme(legend.text=element_text(size=16,face="italic"),legend.title=element_text(size=16))+
  theme(plot.title=element_text(lineheight=1,face="bold",size=19))
ggsave("NMDS_WUnifrac_Proportion_Colour_withClusters_byHIVStatus.pdf")


############## Start Ordination Analysis Multiple Distance measures ##################################
require(plyr)
x.4.0<-xB
#####Intended for de novo OTU-picking only
wh0=genefilter_sample(x.4.0,filterfun_sample(function(x) x>1), A=0.01*nsamples(x.4.0))
x.4.0<-prune_taxa(wh0,x.4.0)
x.4.0 = transform_sample_counts(x.4.0, function(x) ((x/sum(x))))
set.seed(12345)


mymethod="NMDS"
dist=c("bray","canberra","manhattan", "euclidean", "kulczynski", "jaccard", "gower", "altGower", "horn", "mountford" , "binomial")
plist= llply(as.list(dist),function(i,physeq,mymethod){
  set.seed(1)
  ordi=ordinate(physeq,method=mymethod,distance=i,trymax=100)
  plot_ordination(physeq,ordi,"samples",color="RiskGroup")
},x.4.0,mymethod)
names(plist)<-dist

pdataframe = ldply(plist, function(x) {
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "distance"


p = ggplot(pdataframe, aes(Axis_1, Axis_2)) +
  geom_point(size = 2, aes(shape=HIVStatus,color=SexualPractice,fill=SexualPractice)) +
  facet_wrap(~distance, scales = "free")+
  scale_colour_manual(values=c("darkred","darkolivegreen","dodgerblue2")) +
  theme_bw() +
  theme(panel.border=element_blank(),strip.text = element_text(size=12,face="bold"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text.x=element_text(size=4),
        axis.text.y=element_text(size=4))+
  theme(axis.title.x=element_text(size=4),axis.title.y=element_text(size=4))+
  theme(legend.text=element_text(size=16,face="italic"),legend.title=element_text(size=16))+
  stat_ellipse(geom="polygon",aes(fill=SexualPractice),alpha=0.25,level=0.95)+
  ggtitle("Genus level")+
  theme(plot.title=element_text(lineheight=1,face="bold",size=19))
p
ggsave("NMDS_MultipleDistances_WithGenusProportion_WithClusters_byRiskGroup.pdf")
############## End Ordination Analysis Multiple Distance measures ####################################
################ End ordination ########################################


########################################################################
############## Start NMDS cluster Analysis #############################
### Until now we have just mapped our metadata into ordination space and have seen
### that there is an apparently clear association between some of our metadata
### and bacterial composition. Can we further confirm that?
### Let's do a formal clustering analysis

require(vegan)
require(cluster)
x.4.0<-xB
#####Intended for de novo OTU-picking only
wh0=genefilter_sample(x.4.0,filterfun_sample(function(x) x>0), A=0.01*nsamples(x.4.0))
x.4.0<-prune_taxa(wh0,x.4.0)
#x.4.0<-tax_glom(x.4.0,taxrank="Genus")
x.4.0 = transform_sample_counts(x.4.0, function(x) ((x/sum(x))))
set.seed(12345)
x.4.0.ord<-ordinate(x.4.0,"NMDS",distance="bray",trymax=200)
x.4.0.ord

# Distance matrix (Bray-Curtis distance)
Dist<-vegdist(t(otu_table(x.4.0)),method="bray")
# Compute cluster analysis from 2 to 10 groups
################+++++++
Clusterization<-list()
for (i in 1:9){
  Clusterization[[i]]<-pam(Dist,k=i+1)
}

# Silhouette coefficients
Sil<-vector()
for (i in 1:9){Sil[i]<-Clusterization[[i]]$silinfo$avg.width}
# Plot Silhouette coefficient
plot(2:10,Sil,xlab="Number of clusters",ylab="Silhouette coefficient",
     pch=16,lwd=2,type="o",main="Silhouette coefficient")

# So, (by previous analysis) 2 clusters is the best number of groups according to
# silhouette coefficient
Cl<-pam(Dist,k=2)


# Cluster assignation, we define a cluster-label according to PAM analysis
Vector.K2<-as.factor(Cl$clustering)
my.sample_data<-data.frame(sample_data(x.4.0))
my.sample_data$Cluster<-Vector.K2
my.sample_data$SexualPractice<-as.factor(my.sample_data$SexualPractice)
sample_data(xB)<-my.sample_data
data<-data.frame(sample_data(xB))

### Statistics for association between PAM clustering and SexualPractice, in this case
xtabs(~data$Cluster+data$SexualPractice)
fisher.test(xtabs(~data$Cluster+data$SexualPractice))

### Statistics for association between PAM clustering and HIVStatus, in this case
xtabs(~data$Cluster+data$HIV_Status)
fisher.test(xtabs(~data$Cluster+data$HIV_Status))

### Let's use cluster labels for plotting
x.4.0<-xB
wh0=genefilter_sample(x.4.0,filterfun_sample(function(x) x>0), A=0.01*nsamples(x.4.0))
x.4.0<-prune_taxa(wh0,x.4.0)
#x.4.0<-tax_glom(x.4.0,taxrank="Genus")
x.4.0 = transform_sample_counts(x.4.0, function(x) ((x/sum(x))))
set.seed(12345)
x.4.0.ord<-ordinate(x.4.0,"NMDS",distance="bray",trymax=200)

p.4.0.samples=plot_ordination(x.4.0,x.4.0.ord)
p.4.0.samples  + geom_point(size = 5,aes(color=Cluster,shape=SexualPractice)) + ggtitle("Unconstrained PCoA(Bray-Curtis)")+
  scale_colour_manual(values=c("darkgreen","darkorange"))+
  theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),axis.text.x=element_text(size=16),axis.text.y=element_text(size=16))+
  theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
  theme(legend.text=element_text(size=16,face="italic"),legend.title=element_text(size=16))+
  stat_ellipse(geom="polygon",alpha=0.25,aes(fill=Cluster),level=0.95)+
  ggtitle("Genus level")+
  theme(plot.title=element_text(lineheight=1,face="bold",size=19))
ggsave("NMDS_Bray_byNMDSCluster_color.pdf")
capture.output(file="NMDS_Bray_byNMDSCluster_FisherTest.txt",fisher.test(ftable(xtabs(~sample_data(x.4.0)$SexualPractice+sample_data(x.4.0)$Cluster))))

########################################################################
################# END PAM clustering analysis ##########################


############################ Start of piecharts on NMDS #########################################
#################################################################################################
require(mapplots)
x.4.0<-xB
wh0=genefilter_sample(x.4.0,filterfun_sample(function(x) x>1), A=0.01*nsamples(x.4.0))
x.4.0<-prune_taxa(wh0,x.4.0)
x.4.0 = transform_sample_counts(x.4.0, function(x) ((x/sum(x))))
x.4.0<-tax_glom(x.4.0,taxrank="Genus")


#My "Genus" list
my.genus<-c("Prevotella","Bacteroides","Ruminococcus","Catenibacterium")
mytax<-data.frame(tax_table(x.4.0))
my.otus<-rownames(mytax[mytax$Genus %in%my.genus,])
PropData<-data.frame(t(otu_table(x.4.0)))
PropData.2G<-PropData[,my.otus]
PropData.2G<-data.matrix(PropData.2G)

# Minimum and maximum radius
rmin<-0.03; rmax<-0.08
# Sum of 4G
Sum<-apply(PropData[,my.otus],1,sum)

PropData.2G<-PropData.2G/Sum
# Modify Sum (to the interval [rmin,rmax])
Sum2<-(Sum-min(Sum))*(rmax-rmin)/(max(Sum)-min(Sum)) + rmin

# Extract coordinates for a 2-D plot
set.seed(1)
Dist.NMDS<-metaMDS(PropData,k=2,dist="bray")
Puntos<-Dist.NMDS$points


pdf("PiechartsOverNMDS_Bray_BactPrevRumiCate.pdf")
# Plot with the NMDS limits
par(fig=c(0,1,0,.9))
plot(NA,NA,xlim=c(-1.25,1.25),ylim=c(-0.8,1.2),xlab="",ylab="")
# Add text
mtext("NMDS1", side=1, line=3)
mtext("NMDS2", side=4, line=1)

Colours<-c("khaki1","cornflowerblue","brown2","darkgreen")
# Add pies
for (i in 1:nrow(Puntos)){
  #if(data.frame(sample_data(x.4.0))[i,"SexualPractice"]=="HET"){
      add.pie(z=PropData.2G[i,],x=Puntos[i,1],y=Puntos[i,2],radius=Sum2[i],labels=NA,col=Colours)
  #}
}

# Legend (for bar 3)
par(fig=c(0.1,0.9,0.5,1),new=T)
# Empty plot
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend("top",legend=my.genus, fill=Colours,cex=0.8,bty = "n")

dev.off()


############################ Piecharts on NMDS for many genera #################################
################################################################################################

require(mapplots)
x.4.0<-xB
#x.4.0<-subset_samples(x.4.0,SexualPractice=="HET") ### Activate or deactivate accordingly
x.4.0 = transform_sample_counts(x.4.0, function(x) ((x/sum(x))))
x.4.0<-tax_glom(x.4.0,taxrank="Genus")


#My "Genus" list
threshold<-0.02
my.otus<-rownames(mytax[as.vector(rowMeans(otu_table(x.4.0))>threshold),])
my.genus<-paste(mytax[as.vector(rowMeans(otu_table(x.4.0))>threshold),]$Family,"_",mytax[as.vector(rowMeans(otu_table(x.4.0))>threshold),]$Genus,sep="")
PropData<-data.frame(t(otu_table(x.4.0)))
PropData.2G<-PropData[,my.otus]
PropData.2G<-data.matrix(PropData.2G)

# Minimum and maximum radius
rmin<-0.03; rmax<-0.08
# Sum of 4G
Sum<-apply(PropData[,my.otus],1,sum)

PropData.2G<-PropData.2G/Sum
# Modify Sum (to the interval [rmin,rmax])
Sum2<-(Sum-min(Sum))*(rmax-rmin)/(max(Sum)-min(Sum)) + rmin

# Extract coordinates for a 2-D plot
set.seed(1)
Dist.NMDS<-metaMDS(PropData,k=2,dist="bray")
Puntos<-Dist.NMDS$points

pdf("PiechartsOverNMDS_Bray_GeneraOver2.pdf")
par(fig=c(0,3,0,1.2))
plot(NA,NA,xlim=c(-1.25,1.25),ylim=c(-0.8,1.2),xlab="",ylab="")
# Add text
mtext("NMDS1", side=1, line=3)
mtext("NMDS2", side=4, line=1)

require(RColorBrewer)
Colours<-brewer.pal(21,"Spectra")
# Add pies
for (i in 1:nrow(Puntos)){
  if(data.frame(sample_data(x.4.0))[i,"SexualPractice"]=="MSM"){
    add.pie(z=PropData.2G[i,],x=Puntos[i,1],y=Puntos[i,2],radius=Sum2[i],labels=NA,col=Colours)
  }
}

# Legend (for bar 3)
par(fig=c(0.1,0.9,0.5,1),new=T)
# Empty plot
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend("topleft",legend=my.genus, fill=Colours,cex=0.5,bty = "n",ncol=3)
dev.off()
########################### Stop of piecharts on NMDS ##########################################
################################################################################################



#######Adonis TESTs##########################################
### Until now we have addressed association

x.6.0<-xB
wh0=genefilter_sample(x.6.0,filterfun_sample(function(x) x>1), A=0.01*nsamples(x.6.0))
x.6.0<-prune_taxa(wh0,x.6.0)
x.6.0.genus<-tax_glom(x.6.0,taxrank="Genus")
x.6.0<-x.6.0.genus
x.6.0 = transform_sample_counts(x.6.0, function(x) (((1* (x))/sum(x))))
#### Non-parametric analysis of variance using NPMANOVA through adonis function {vegan}package

x.6.0.response<-t(otu_table(x.6.0))
metadata<-data.frame(sample_data(x.6.0))

x.6.0.explanatory<-metadata[,c("SexualPractice","RiskGroup","HIVStatus","Age")]
for(var in colnames(x.6.0.explanatory)){
  explanatory=data.frame(x.6.0.explanatory[,eval(var)])
  myadonis<-adonis(x.6.0.response~.,data=explanatory)
  capture.output(paste("Adonis on: ",var),file="adonis_tests.txt",append=T)
  capture.output(myadonis,file="adonis_tests.txt",append=T)
}

adonis(x.6.0.response~.,data=x.6.0.explanatory)
x.6.0.adonis<-adonis(x.6.0.response~.,data=x.6.0.explanatory)
simper(x.6.0.response,x.6.0.explanatory[,"HIVStatus"])
simper(x.6.0.response,x.6.0.explanatory[,"SexualPractice"])

######################################### END of Adonis Tests #######################################
#####################################################################################################


#### Heatmap and Dendrograms
x.4.1<-xB
x.4.1<-tax_glom(x.4.1,taxrank="Genus")
x.4.1 = transform_sample_counts(x.4.1, function(x) ((x/sum(x))))

# A function to draw a heatmap containg bacterial genus abundance proportions
myHeatmap<-function(xx,level="Genus",distance="wunifrac",filter=NULL){
  require("phyloseq")
  require("RColorBrewer")
  scaleyellowred <- colorRampPalette(c( "steelblue","darkred"), space = "rgb")(100)
  data.dist<-UniFrac(xx,weighted=T,normalized=T)
  col.clust<-hclust(data.dist,"aver")
  xx<-prune_taxa(filter,xx)
  provat<-data.frame(tax_table(xx)[colnames(as.matrix(t(otu_table(xx)))),])
  myOtuTable<-as.matrix(t(otu_table(xx)))
  TaxFields<-c("Order","Family","Genus")
  provat$TaxID<-do.call("paste",c(provat[TaxFields],sep="_"))
  colnames(myOtuTable)<-provat$TaxID
  #par?par(paper="a4r")
  heatmap((myOtuTable),Rowv=as.dendrogram(col.clust),col=scaleyellowred,margins = c(13,13),size=1)
}

#Create Filter to select taxa represented in heatmap
# At least 1% presence in at least two samples

wh4.1=genefilter_sample(x.4.1,filterfun_sample(function(x) x>0.01), A=2)
pdf("Heatmap_WUnifracDendrogam_GenusGlom_1perc_2.pdf",paper="A4r")
myHeatmap(x.4.1,filter=wh4.1)
dev.off()

myPheatmap<-function(xx,level="Genus",distance="wunifrac",filter=NULL){
  require("phyloseq")
  require("RColorBrewer")
  require("pheatmap")
  require("gplots")
  scaleyellowred <- colorRampPalette(c( "steelblue","darkred"), space = "rgb")(100)
  data.dist<-UniFrac(xx,weighted=T,normalized=F)
  col.clust<-hclust(data.dist,"aver")
  plot(col.clust)
  xx<-prune_taxa(filter,xx)
  provat<-data.frame(tax_table(xx)[colnames(as.matrix(t(otu_table(xx)))),])
  myOtuTable<-as.matrix(t(otu_table(xx)))
  TaxFields<-c("Order","Family","Genus")
  provat$TaxID<-do.call("paste",c(provat[TaxFields],sep="_"))
  colnames(myOtuTable)<-provat$TaxID
  #par?par(paper="a4r")
  pheatmap((myOtuTable),clust_cols=T,clust_rows=T,clustering_distance_rows=data.dist,clustering_method="average")
}
pdf("PHeatmap_WUnifracDendrogam_GenusGlom_1perc_2.pdf",paper="A4r")
myPheatmap(x.4.1,filter=wh4.1)
dev.off()

####
#### Numeric Boxplots
compareTwoGroups<-function(mydata=NULL,variable=NULL,category1=NULL,category2=NULL,fileForPlot=NULL,minCounts=500,maxAlpha=0.01, design=NULL){
  fileForPlot=paste("NegBin_DiffTest_",as.character(deparse(substitute(mydata))),"_",variable,"_",category1,"vs",category2,".pdf",sep="")
  print(paste("Output file for Plots: ",fileForPlot))
  fileForTable=paste("NegBin_DiffTest_",as.character(deparse(substitute(mydata))),"_",eval(variable),"_",eval(category1),"vs",eval(category2),".txt",sep="")
  print(paste("Output file for Table: ",fileForTable))
  stringForTitle=paste(variable,"/",category1,"vs",category2)
  kostic <- mydata
  myTaxTable<-tax_table(mydata)
  LastRank<-rank_names(mydata)[length(rank_names(mydata))][[1]]
  myLastRankName<-sort(unique(myTaxTable[,"Phylum"]))
  print(myLastRankName)
  require(phyloseq)
  require(DESeq2)
  kostic <- prune_samples(sample_sums(kostic) > minCounts, kostic)
  diagdds = phyloseq_to_deseq2(kostic, as.formula(design))
  #print(diagdds)
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(diagdds), 1, gm_mean)
  diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
  #colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("untreated","treated"))
  #colData(diagdds)$condition<-factor(colData(diagdds)$condition,levels=c(category1,category2))
  diagdds = DESeq(diagdds, fitType="parametric",test="Wald")
  res=results(diagdds,contrast=c(variable,category1,category2))
  #res=results(diagdds)
  print(res)
  res = res[order(res$padj, na.last=NA), ]
  sigtab = res[(res$padj < 1), ]
  #   sigtab = res[(res$padj < maxAlpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(kostic)[rownames(sigtab), ], "matrix"))
  sigtab$OtuID<-rownames(sigtab)
  #   head(sigtab)
  write.table(sigtab,file=fileForTable,sep="\t")
  sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=levels(sort(sigtab$Phylum)))
  #   sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=levels(myLastRankName))
  print(dim(sigtab))
  #Cleanup for Positive enrichment in csigtabarcinoma
  posigtab=sigtab
  #posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
  #posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
  posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", rank_names(mydata))]
  library("ggplot2")
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set3", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  sigtabgen=sigtab
  sigtabgen = subset(sigtab, !is.na(Genus))
  sigtabgen = subset(sigtabgen, padj <= maxAlpha)
  sigtabgen = subset(sigtabgen, baseMean >= 10)
  sigtabgen = subset(sigtabgen, sigtabgen$Genus != "unclassified")
  # Phylum order
  #   x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
  #   x = sort(x, TRUE)
  #   sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
  #   sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=sort(as.character(sigtabgen$Phylum)))

  # Genus order
  if(as.character(LastRank)== "Genus"){
    sigtabgen$LastRank<-sigtabgen[,"Genus"]
  }else{
    sigtabgen$LastRank<-paste(sigtabgen[,"Genus"]," ",sigtabgen[,as.character(LastRank)])
  }
  x = tapply(sigtabgen$log2FoldChange, sigtabgen$LastRank, function(x) min(x))
  x = sort(x, TRUE)

  sigtabgen$LastRank = factor(as.character(sigtabgen$LastRank), levels=names(x))
  sigtabgen$log2Counts<-sigtabgen$baseMean
  sigtabgen$alpha<- 1 - sigtabgen$padj
  #pdf(fileForPlot)
  p<-ggplot(sigtabgen,aes(x=LastRank,y=log2FoldChange))
  p<-p+geom_point(aes(colour=Phylum,size=baseMean,alpha=alpha))
  p<-p+ guides(colour = guide_legend(override.aes = list(size=10)))
  p<-p+scale_size_continuous(range=c(10,20))
  #+geom_point(aes(size=sigtabgen$log2Counts))+scale_size_area()
  p<-p+theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=0.5,size=10))
  #p<-p+theme(legend.key.size=unit(1,"cm"))
  p<-p+ ggtitle(paste(stringForTitle," Data:",as.character(deparse(substitute(mydata))))) +
    theme(plot.title = element_text(lineheight=.7, face="bold"))+coord_flip()+
    theme(axis.text.y = element_text( size=16)) + geom_hline(xintercept=0,colour="darkred", linetype = "longdash")
  print(p)
  #ggplot(sigtabgen, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6)+scale_size(range=c(1,5))+
  # theme(axis.text.x = element_text(angle = -90, hjust = 0,size=3, vjust=0.5), legend.key.size=unit(0.5,"cm"),legend.text=element_text(size=3))
  ggsave(filename = fileForPlot,dpi=600,width=11, height=8.5)
  return(sigtab)
  #dev.off()
}
statisticValue<-function(value=NULL){
  if(value>0.1){
    return (paste("=",round(value,digits=2)))
  }
  else if(value>=0.05){
    return (paste("=",round(value,digits=2)))
  }
  else if(value>=0.01){
    return (paste("<0.05"))
  }
  else if(value>=0.001){
    return (paste("<0.01"))
  }
  else if(value>=0.0001){
    return(paste("<0.001"))
  }
  else{return (paste("<0.0001"))}
}

x.3.0<-xB
compareTwoGroups(mydata=x.3.0,variable="HIVStatus",category1="HIVpos",category2="HIVneg",design=~HIVStatus,maxAlpha=0.1)

#################### Additional section


### In the following section we will just output the data in a specific format
### for additional analysis using external tools. This will create some files that
### can be imported downstres
wh0=genefilter_sample(xB,filterfun_sample(function(x) x>1), A=0.01*nsamples(xB))
xBF<-prune_taxa(wh0,xB)
####Format Output for LEFSe. We use some pre-coded function for this.
myLefseFunction<-function(psObject,rank){
  otus<-data.frame(otu_table(tax_glom(psObject,taxrank=rank)))
  taxons<-data.frame(tax_table(tax_glom(psObject,taxrank=rank)))
  for (j in 1:dim(taxons)[1]){
    taxid<-taxons[j,1]
    for (i in 2:dim(taxons)[2]){
      if(!(is.na(taxons[j,i]))){
        taxid<-paste(taxid,taxons[j,i],sep="|")
      }
    }
    taxons[j,"Taxon"]<-taxid
  }
  taxons$OtuID<-rownames(taxons)
  taxons<-(taxons[,c("OtuID","Taxon")])
  otus$OtuID<-rownames(otus)
  taxons<-merge(taxons,otus,by="OtuID")
  taxons$OtuID<-NULL
  return(taxons)
}
myPhyloseqToLefse<-function(psObject){
  xx<- transform_sample_counts(psObject, function(x) ((( x)/sum(x))))
  xx.genus<-myLefseFunction(xx,rank="Genus")
  xx.family<-myLefseFunction(xx,rank="Family")
  xx.order<-myLefseFunction(xx,rank="Order")
  xx.class<-myLefseFunction(xx,rank="Class")
  xx.phylum<-myLefseFunction(xx,rank="Phylum")
  xx.kingdom<-myLefseFunction(xx,rank="Kingdom")
  xx<-rbind(xx.genus,xx.family,xx.order,xx.class,xx.phylum,xx.kingdom)
  rownames(xx)<-xx$Taxon
  xx$Taxon<-NULL
  xx<-t(xx)
  rownames(xx)<-gsub("X","",rownames(xx))
  rownames(xx)<-gsub("\\.","-",rownames(xx))
  xx<-merge(data.frame(sample_data(psObject)),xx,by="row.names")
  xx$Row.names<-NULL
  xx<-t(xx)
  colnames(xx)<-xx[1,]
  xx<-data.frame(xx[-1,])
  xx$SampleID<-rownames(xx)
  xx<-xx[,c(ncol(xx),1:(ncol(xx)-1))]
  return(xx)
}
xx<-myPhyloseqToLefse(xBF)
write.table(xx,paste(as.character(date()),"Lefse.txt"),sep="\t",row.names=F)
### End LEFSE #############################

### Format Output for STAMP #######
myStampFunction<-function(psObject){
  otu<-data.frame(otu_table(psObject))
  tax<-data.frame(tax_table(psObject))
  tax<-data.frame(merge(tax,otu,by="row.names"))
  tax$Row.names<-NULL
  return(tax)
}
xx<-myStampFunction(xBF)
write.table(xx,paste(as.character(date()),"STAMP.txt"),sep="\t",row.names=F)
### End format for STAMP ####



