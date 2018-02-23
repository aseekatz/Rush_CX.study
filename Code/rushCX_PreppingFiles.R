#### Rush CX study: Code used for manuscript
##### Anna M. Seekatz
##### 5.20.17


###### Raw data file used in analyses:
	- sample information and de-identified metadata: rushCX.meta_strict.CRE.txt
	- analytic sample data: SMB1AnalyticDataset.txt
	- OTU counts (filtered in mothur): rushfinal.CXsamples.0.03.0.03.filter.shared
	- OTU classifications: mothurfiles_rushfinal/rushfinal.final.0.03.cons.taxonomy
		- note: this includes all OTUs processed together with a larger study 
	- summary files calculated in mothur:
		- rushfinal.CXsamples.0.03.thetayc.0.03.lt.nmds.axes
		- rushfinal.CXsamples.0.03.thetayc.0.03.lt.pcoa.axes
		- rushfinal.CXsamples.0.03.groups.summary
	- phylotype file (up to genus-level): rushfinalCX_wang.tax.summary.txt
		- note: this included counts for a larger data set, and has been filtered
	- pairwise calculations of different beta diversity measures: rushfinal.CXsamples.0.03.summary
		- this was combined with metadata 
		
###### Data files merged with meta data or manipulated otherwise:
	- rushCX_filtered.shared.frac.txt --> relative abundance of OTUs
	- rushfinalCX_summary.txt --> metadata with summary scores
		- rushfinalCX_summary.w.pam.txt: additional summary file with pam clustering
	- rushCX.03_taxnames.txt --> specific OTU classifications
	- rushfinalCX_genera2pfrac.txt, rushfinalCX_all.genera.txt --> phylotype files 
	- rushfinal.CXsamples.0.03.summary, rushfinalCX_summary.w.pam.txt --> beta diversity pairwise calcuations, with relevant metadata

####### OTU relative abundance:

```
# read in files and convert to relative abundance:
shared<-read.table(file="mothurfiles/rushfinal.CXsamples.0.03.0.03.filter.shared", header=TRUE, row.names=2)
cx<-read.table(file="Data/rushfinalCX_summary.txt", header=TRUE, sep="\t")
names<-as.character(cx$group)		# these are all of the subset samples
shared.good<-shared[which(rownames(shared) %in% names), 3:ncol(shared)]
dd.all<-as.matrix(shared.good)
sub.dd.frac<-dd.all/rowSums(dd.all)
#write.table(sub.dd.frac, file="Data/rushCX_filtered.shared.frac.txt", sep="\t", quote=FALSE, col.names=NA)

# note: this FILTERED file, although it has many OTUs, is very selective...in some samples, up to 85% of the sequences have been removed
# let's look at the full file a bit more closely:
rawshared<-read.table(file="mothurfiles/rushfinal.CXsamples.0.03.shared", header=TRUE, row.names=2)
cx<-read.table(file="Data/rushfinalCX_summary.txt", header=TRUE, sep="\t")
names<-as.character(cx$group)		# these are all of the subset samples
shared.good<-rawshared[which(rownames(rawshared) %in% names), 3:ncol(rawshared)]
shared.good<-shared.good[rowSums(shared.good) > 0,]

# get rid of all the samples 
min(rowSums(shared.good))	#10421 is the minimum--pretty good

# let's apply filtering within R
library(vegan)
set.seed(1984)	# to get consistent results
sub_size <- as.numeric(min(rowSums(shared.good)))
	# Rarefy shared file to appropriate min (subsample):
shared <- as.data.frame(t(shared.good))

for (index in 1:ncol(shared)){
  shared[,index] <- t(rrarefy(shared[,index], sample=sub_size))}
rm(index, sub_size)

# remove OTUs that are no longer existent due to rarification (0 totals):
#colSums(shared)		# if you want to check that everything was rarefied
shared<-shared[rowSums(shared) > 0,]
dim(shared)
	# this reduced OTUs to 6540

# last step, if you want, is to filter by prevalence
# Filter out columns that have values in at least 2 samples (ignores first column if needed)
filter_table <- function(data) {
  drop <- c()
  if (class(data[,1]) != 'character') {
    if (sum(data[,1] != 0) < 2) {
      drop <- c(drop, colnames(data)[1])
    }
  }
  for (index in 2:ncol(data)) {
    if (sum(data[,index] != 0) < 2) {
      drop <- c(drop, colnames(data)[index])
    }
  }
  filtered_data <- data[,!(colnames(data) %in% drop)]
  filtered_data$rareOTUs <- rowSums(data[,(colnames(data) %in% drop)])
  return(filtered_data)
}
filtered.shared <- as.data.frame(t(shared))
filtered.shared <- filter_table(filtered.shared) # OTU set reduced to 561--this is probably more appropriate than the method used in mothur, since all the samples are so different from each other
dim(filtered.shared)
rowSums(filtered.shared)
	# using at least 2 samples reduced the data set to 1182
	# using at least 3 samples reduced the data set to 832
	# let's look at the weird OTUs in some of these...
filtered.otus<-colnames(filtered.shared)
shared2<-as.data.frame(t(shared))
test<-shared2[,!(colnames(shared2) %in% filtered.otus)]
test$remainder<-rowSums(shared2[,(colnames(shared2) %in% filtered.otus)])
relabund<-as.matrix(test/rowSums(test))
relabund<-as.data.frame(relabund)
relabund<-relabund[order(relabund$remainder),]
relabund<-subset(relabund, select=-c(remainder))
max.relabund<-relabund[,colSums(relabund) > 0.01]
weirdOTUs<-as.character(colnames(max.relabund))
	# read in taxonomic classifications:
tax<-read.table(file="manuscript/Git/Data/rushfinal.03_taxnames_allOTUs.txt", header=TRUE)
weirdtax<-tax[tax$OTU %in% weirdOTUs, ]
	# all of these look like they are likely not contamination (although we do not know)
	# since these are all cross-sectional, not sure what to do about weird OTUs that only show up in 1 patient...
	# we have this file if needed--for now, we can still use the mothur file
	# will output this just in case it is needed

# convert to relative abundance:
otufrac<-as.data.frame(as.matrix(filtered.shared)/rowSums(as.matrix(filtered.shared))*100)
#write.table(otufrac.meta, file="rushCX_otufrac.moreOTUs.txt", sep="\t", quote=FALSE, col.names=NA)

```

####### summary file with metadata:

```
# read in files and merge:

# since this was taken from a large group file, let's redo the analysis in mothur and get a smaller file:
meta<-read.csv("Data/rushCX.meta_strict.CRE.txt", header = TRUE, sep = "\t", strip.white=TRUE, na.strings = c("", " "))
# 0.03:
pcoa3<-read.table(file="mothurfiles/rushfinal.CXsamples.0.03.thetayc.0.03.lt.pcoa.axes", header=TRUE)
	pcoa3<-pcoa3[,1:4]
	colnames(pcoa3)[2:4] <- paste("pcoa03", colnames(pcoa3)[2:4], sep = "_")
	colnames(pcoa3)[1]<-"sampleID"
nmds3<-read.table(file="mothurfiles/rushfinal.CXsamples.0.03.thetayc.0.03.lt.nmds.axes", header=TRUE)
	nmds3<-nmds3[1:4]
	colnames(nmds3)[2:4] <- paste("nmds03", colnames(nmds3)[2:4], sep = "_")
	colnames(nmds3)[1]<-"sampleID"
sum3<-read.table(file="mothurfiles/rushfinal.CXsamples.0.03.groups.summary", header=TRUE)
	sum3<-subset(sum3, select=-c(label))
	colnames(sum3)[2:16] <- paste(colnames(sum3)[2:16], "03", sep = "_")
	colnames(sum3)[1]<-"sampleID"

combined.pcoa<-merge(meta, pcoa3, by.x=c("group"), by.y=c("sampleID"), all.x=TRUE)
combined.nmds<-merge(combined.pcoa, nmds3, by.x=c("group"), by.y=c("sampleID"))
combined.sum<-merge(combined.nmds, sum3, by.x=c("group"), by.y=c("sampleID"))
#write.table(combined.sum, 'Data/rushfinalCX_summary.txt',quote=FALSE,sep="\t", col.names=TRUE, row.names=FALSE)


```

####### Taxonomy file from mothur:

```

# 0.03 taxonomy file:
taxonomy_file<-read.table(file="mothurfiles/rushfinal.final.0.03.cons.taxonomy", header=TRUE)
tax <- taxonomy_file$Taxonomy
tax <- gsub("\\(\\d*\\)", "", tax)
tax <- gsub(";unclassified", "", tax)
tax <- gsub("_1", "", tax)
tax <- gsub(";$", "", tax)
tax <- gsub("/.*", "", tax)
tax <- gsub(".*;", "", tax)
tax.names <-paste(taxonomy_file$OTU, tax, sep="_")
#tax.names <-gsub("000", "", tax.names)
taxonomy_file$taxname<-tax.names
phylum <- taxonomy_file$Taxonomy
phylum <- gsub("\\(\\d*\\)", "", phylum)
phylum <- gsub("Bacteria;", "", phylum)
phylum <- gsub(";$", "", phylum)
phylum <- gsub(";.*", "", phylum)
taxonomy_file$phylum<-phylum
tax<-taxonomy_file
tax<-tax[order(tax$OTU),]
#write.table(tax, file="Data/rushfinal.03_taxnames_allOTUs.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# then, filter to get OTUs identified in specific CX study set (after filtering)
shared<-read.table(file="mothurfiles/rushfinal.CXsamples.0.03.0.03.filter.shared", header=TRUE)
tax03<-read.table(file="mothurfiles/rushfinal.0.03.taxonomy.names.txt", header=TRUE)
OTUs<-as.character(colnames(shared)[4:ncol(shared)])
filtered.tax03<-tax03[tax03$OTU %in% OTUs, ]
#write.table(filtered.tax03, file="Data/rushCX.03_taxnames.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

```

###### Phylotype file in mothur:

```
# read in the taxonomy file, and filter to rushCX samples:
taxo<-read.table(file="mothurfiles/new_taxonomy/rushfinal.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.tax.summary", header=TRUE)
taxo.samples<-taxo[, 6:ncol(taxo)]
cx<-read.table(file="Data/rushfinalCX_summary.txt", header=TRUE, sep="\t")
names<-as.character(cx$group)
taxo.names<-taxo.samples[, colnames(taxo.samples) %in% names]	#get rushCX samples only
taxo.filtered<-cbind(taxo[1:5], taxo.names)
taxo.present<-taxo.filtered[rowSums(taxo.filtered[6:ncol(taxo.filtered)]) > 0, ]
#write.table(taxo.present, file="mothurfiles/rushfinalCX_wang.tax.summary.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)	

# filter out genus-level assignments only, and assign a phylum:
tax<-taxo.present
tax2<-tax[which(tax$taxlevel==2), ]
tax2[, c("rankID", "taxon")]
tax6<-tax[which(tax$taxlevel==6), ]
tax6$rankID<-gsub("^0.1.2.*", "20_Euryarchaeota", tax6$rankID)
tax6$rankID<-gsub("^0.1.3.*", "20_Thaumarchaeota", tax6$rankID)
tax6$rankID<-gsub("^0.2.1\\..*", "10_Acidobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.2\\..*", "04_Actinobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.4\\..*", "11_Bacteria_unclassified", tax6$rankID)
tax6$rankID<-gsub("^0.2.5\\..*", "01_Bacteroidetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.6\\..*", "20_Candidatus_Saccharibacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.7\\..*", "20_Chlamydiae", tax6$rankID)
tax6$rankID<-gsub("^0.2.9\\..*", "20_Cyanobacteria/Chloroplast", tax6$rankID)
tax6$rankID<-gsub("^0.2.8\\..*", "20_Chloroflexi", tax6$rankID)
tax6$rankID<-gsub("^0.2.10..*", "20_Deferribacteres", tax6$rankID)
tax6$rankID<-gsub("^0.2.13..*", "02_Firmicutes", tax6$rankID)
tax6$rankID<-gsub("^0.2.14..*", "06_Fusobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.15..*", "20_Gemmatimonadetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.17..*", "20_Ignavibacteriae", tax6$rankID)
tax6$rankID<-gsub("^0.2.18..*", "20_Lentisphaerae", tax6$rankID)
tax6$rankID<-gsub("^0.2.20..*", "20_Nitrospirae", tax6$rankID)
tax6$rankID<-gsub("^0.2.21..*", "20_Parcubacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.22..*", "20_Planctomycetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.23..*", "03_Proteobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.24..*", "09_Spirochaetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.25..*", "08_Synergistetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.26..*", "07_Tenericutes", tax6$rankID)
tax6$rankID<-gsub("^0.2.27..*", "05_Verrucomicrobia", tax6$rankID)
tax6$rankID<-gsub("^0.3.1..*", "11_unknown_unclassified", tax6$rankID)
colnames(tax6)[2]<-"phylum"

# filter some columns, turn into a matrix and check for duplication:
subtax6<-subset(tax6, select=-c(taxlevel, daughterlevels))
subtax6<-subtax6[order(subtax6$phylum, -subtax6$total), ]
taxmatrix<-subtax6[, c(4:ncol(subtax6))]
which(duplicated(subtax6$taxon))
	# fix the duplicated row:
subtax6$taxon<-as.character(subtax6$taxon)
subtax6$taxon[268]<-"Actinobacteria_unclassified2"
subtax6$taxon<-as.factor(subtax6$taxon)
rownames(taxmatrix)<-subtax6$taxon
genera<- taxmatrix[, colSums(taxmatrix)>5000,]		# should already be filtered
	# get rel. abund fraction:
genmatrix<-as.data.frame(t(genera))
genera.fr<-genmatrix/rowSums(genmatrix)*100
genus.fr<-t(genera.fr)
all.genera<-cbind(subtax6[1:3], genus.fr)
#write.table(all.genera, file="Data/rushfinalCX_all.genera.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
	# get top 1%:
phyla<-subtax6[1:3]
genus1<- genus.fr[rowSums(genus.fr>=1)>=1,]
namelist<-as.character(rownames(genus1))
phyla1p<-phyla[phyla$taxon %in% namelist, ]
genera1<-cbind(phyla1p, genus1)
	# get top 2%
genus2<- genus.fr[rowSums(genus.fr>=2)>=2,]
namelist<-as.character(rownames(genus2))
phyla2p<-phyla[phyla$taxon %in% namelist, ]
genera2<-cbind(phyla2p, genus2)
#write.table(genera2, file="Data/rushfinalCX_genera2pfrac.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# note: later, pam clustering was added to this to produce rushfinalCX_summary.w.pam.txt (see Fig3 code)


```

###### Pairwise beta diversity measurements, with metadata:

```

# read in files and add relevant meta info:
# (done in CX_PreppingFiles)

meta<-read.table(file="/Data/rushfinalCX_summary.w.pam.txt", header=TRUE, sep="\t")
mmeta<-meta[, c("group", "SAMStudy_Number", "sample_set", "cre_status", "phylopam3_cluster")]
comp <- read.csv("mothurfiles/rushfinal.CXsamples.0.03.summary",  stringsAsFactors=FALSE, header=T, sep="\t") 
#comp <- subset(comp, select=-c(X, label))

# since this file is so large, let's filter out some data:
mmeta<-meta[, c("group", "SAMStudy_Number", "sample_set", "cre_status", "phylopam3_cluster", "Antibiotic_all")]
data.names<-as.character(mmeta$group)
compf<-comp[comp$s1 %in% data.names & comp$s2 %in% data.names, ]

# read in file and merge with relevant meta
m<-merge(mmeta, compf, by.x=c("group"), by.y=c("s1"))
s1names<-paste("s1", colnames(m[1:6]), sep = "_")
allnames<-c(s1names, colnames(m[7:15]))			
colnames(m) <- allnames
m2<-merge(mmeta, m, by.x=c("group"), by.y=c("s2")) 			
s2names<-paste("s2", colnames(m2[1:6]), sep = "_")
allnames<-c(s2names, s1names, colnames(m[8:15]))			
colnames(m2) <- allnames
#write.table(m2, file="Data/CX_alldist.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)	


```