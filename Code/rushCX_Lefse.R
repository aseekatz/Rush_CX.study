## 091516 Rush CX analysis

##Using unfiltered shared for amova
setwd ('/Users/christinebassis/Desktop/091516_Rush_CX/')
shared<-read.table(file='rushfinal.CXsamples.0.03.shared', header=TRUE, row.names=2) 
dim(shared)
meta<-read.table("rushCX_meta.txt", header = TRUE, sep = "\t", strip.white=TRUE, na.strings = c("", " "))
dim(meta)
rushCX.shared.wmeta<-merge(shared, meta, by.x=c("row.names"), by.y=c("seqID"))
write.table(rushCX.shared.wmeta, "rushCX.shared.wmeta.txt", quote=FALSE,sep="\t", col.names=NA)

##opened rushCX.shared.wmeta.txt in Excel and moved metadata columns to the front and saved as rushCX.shared.wmeta_first.txt
  #sorted by cre_status, then by cre_org_category
  #kept only if: cre_status=KPC POSITIVE and cre_org_category=KLEB (n=32)
          #  or cre_status=negative (n=99)
          # do we care how long they have been in LTACH?
#saved as rushCX.shared.wmeta_first_crepos_kleb_or_neg_only.txt
  #deleted metadata columns, labeled "Group" column and moved to 2nd column position so it goes: label, Group, numOTUs, Otu000001.....
  #saved as rushCX.shared_crepos_kleb_or_neg_only.txt
#made design file with cre_status (changed KPC POSITIVE to KPC_POSITIVE) saved as cre_status_strict.design.txt

#in mothur v.1.38.1 (Mac version)
#mothur > dist.shared(shared=/Users/christinebassis/Desktop/091516_Rush_CX/rushX.shared_crepos_kleb_or_neg_only.txt, calc=thetayc)

#Using 1 processors.
#0.03

#Output File Names: 
#  /Users/christinebassis/Desktop/091516_Rush_CX/rushCX.shared_crepos_kleb_or_neg_only.thetayc.0.03.lt.dist


#mothur > amova(phylip=/Users/christinebassis/Desktop/091516_Rush_CX/rushCX.shared_crepos_kleb_or_neg_only.thetayc.0.03.lt.dist, design=/Users/christinebassis/Desktop/091516_Rush_CX/cre_status_strict.design.txt)
#KPC_POSITIVE-negative	Among	Within	Total
#SS	0.714127	56.6587	57.3728
#df	1	129	130
#MS	0.714127	0.439215

#Fs:	1.62592
#p-value: 0.017*
  
#  Experiment-wise error rate: 0.05
#If you have borderline P-values, you should try increasing the number of iterations

#Output File Names: 
#  /Users/christinebassis/Desktop/091516_Rush_CX/rushCX.shared_crepos_kleb_or_neg_only.thetayc.0.03.lt.amova


###Using filtered for LEfSe
setwd ('/Users/christinebassis1/Box Sync/091516_Rush_CX/')
shared<-read.table(file='Rushfinal.CXsamples.0.03.0.03.filter.shared', header=TRUE, row.names=2) 
dim(shared)
meta<-read.table("rushCX_meta.txt", header = TRUE, sep = "\t", strip.white=TRUE, na.strings = c("", " "))
dim(meta)
rushCX.shared.filtered.wmeta<-merge(shared, meta, by.x=c("row.names"), by.y=c("seqID"))
write.table(rushCX.shared.filtered.wmeta, "rushCX.shared.filtered.wmeta.txt", quote=FALSE,sep="\t", col.names=NA)

mothur > lefse(shared=/Users/christinebassis1/Desktop/091516_Rush_CX/rushCX.shred.filtered.crepos_kleb_or_neg_only.txt, design=/Users/christinebassis1/Desktop/091516_Rush_CX/cre_status_strict.design.txt)

You did not provide a class, using value0.

0.03

Number of significantly discriminative features: 12 ( 51 ) before internal wilcoxon.
Number of discriminative features with abs LDA score > 2 : 12.

Output File Names: 
  /Users/christinebassis1/Desktop/091516_Rush_CX/rushCX.shared.filtered.crepos_kleb_or_neg_only.0.03.lefse_summary