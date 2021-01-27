#######
## Master in Omics Data Analysis, 2018-19 Universitat de Vic
## Metagenomics
## Marc Noguera-Julian, PhD. February, 3rd,4th and 5th, 2020
######

## You can access analysis code at https://github.com/mnoguera/MasterDadesOmiques 
## or as a zip file in the campus. Use GitHub for latest version.

## DataSet for practical drylab work on statistical analysis16s metagenomics data
## We will be analyzing the data set from 110 sample in non-infected and HIV-infected patients
## From the publication:
## Noguera-Julian, M., Rocafort, M., Guillén, Y., Rivera, J., Casadellà, M., Nowak, P., Hildebrand, F., Zeller, G., Parera, M., Bellido, R., et al. (2016). Gut Microbiota Linked to Sexual Preference and HIV Infection. EBioMedicine 5, 135–146.

## We will be analyzing the OTU table from month 0, along with a subset of the available metadata
## You will need the Following files:
##   -  Metadata_Workshop_v1.0.csv: Contains the Metadata
##   -  final.an.unique_list.0.03.cons.taxonomy: Contains the taxonomical annotation of the OTUs in the OTU table
##   -  SampleRenamed.shared: Contains the count table of the OTUs for each sample
##   -  OtuRepSeq.tree: Contains a phylogenetic tree based on the representative sequences for each OTU


## For the analysis we will be using R and the following packages
> packageVersion("phyloseq")  
[1] ‘1.26.1’
> packageVersion("vegan")
[1] ‘2.5.4’
> packageVersion("ggplot2")
[1] ‘3.1.0’
> packageVersion("DESeq2")
[1] ‘1.22.2’
> packageVersion("RColorBrewer")
[1] ‘1.1.2’
> packageVersion("doBy")
[1] ‘4.6.2’
> packageVersion("mapplots")
[1] ‘1.5.1’
