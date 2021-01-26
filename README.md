#######
## Master in Omics Data Analysis, 2020-2021 Universitat de Vic
## Metagenomics
## Marc Noguera-Julian, PhD. Jan, 26th, 27th, 28th 2021
######

DataSet for practical drylab work on 16s metagenomics data
RawData: contains Quality controlled 10k illumina fastq 
			https://drive.google.com/file/d/13OdeTDwXCGjLfivG4ZWpPb_R4D8kiF3z/view?usp=sharing
Config: Contains metadata for further analysis 

We will use SilvaDB to classify our sequence dataset
SilvaDB is the most complete database of rRNA sequences and 
provides reduced and extended alignments for analysis
Check it out at: http://www.arb-silva.de

For the practical session we will use these reference datasets formatted for DADA2
These datasets are quite large and may take quite a long time to download 
depending on your internet connection speed.

Please download all DB from
		https://drive.google.com/file/d/1Js4IVVsgZr6N6UEkm8xW23RNL5tXbGXs/view?usp=sharing	
And decompress into Training Subdir

We will use RDP algorithm by Wang et al (2012)
for sequence classification as implemented in DADA2
Find the publication here http://aem.asm.org/content/73/16/5261

Drylab work will use DADA2 https://benjjneb.github.io/dada2/

Please install within your Rstudio environment:

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2", version = "3.10")

Additional R-packages: phyloseq, shortRead, digest

