##### Metagenomics course, 2015-16
##### Marc Noguera-Julian, PhD, IrsiCaixa, UVic AIDS Chair
##### 2016, Feb 2nd.

require(phyloseq)
##### A Functions file
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
  xx<- transform_sample_counts(x, function(x) ((( x)/sum(x))))
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
  xx<-merge(data.frame(sample_data(x)),xx,by="row.names")
  xx$Row.names<-NULL
  xx<-t(xx)
  colnames(xx)<-xx[1,]
  xx<-data.frame(xx[-1,])
  xx$SampleID<-rownames(xx)
  xx<-xx[,c(ncol(xx),1:(ncol(xx)-1))]
  return(xx)
}