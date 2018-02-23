#### Rush CX study: Fig. 4, relative abundance of differentially abundant OTUs as identified by Lefse
##### Anna M. Seekatz
##### 12.18.17

###### Figures created:
	- Fig. 4: relative abundance of differentially abundant OTUs using lefse across calibration, validation, and all samples


####### Files used:
	- rushCX.shared.filtered.crepos_kleb_or_neg_only.0.03.lefse_summary.w.classifications_all_cal_val_comparison
	- rushCX_filtered.shared.frac.txt
	- rushfinalCX_summary.w.pam.txt

######## Fig. 4: Relative abundance of Lefse-identified OTUs in valibration, validation cohorts

```
library(gplots)


# step 1: create OTU labels
# read in data, and plot with bargraphs similar to script used for RF plotting
lf<-read.table(file="LEfSe/rushCX.shared.filtered.crepos_kleb_or_neg_only.0.03.lefse_summary.w.classifications_all_cal_val_comparison.txt", header=TRUE, sep = "\t", strip.white=TRUE, na.strings = c("", " "))
	# make list of OTUs (by group identified)
all<-as.character(lf[lf$All_Class %in% c("KPC_POSITIVE", "negative"), c("OTU")])
cal<-as.character(lf[lf$Cal_Class %in% c("KPC_POSITIVE", "negative"), c("OTU")])
val<-as.character(lf[lf$Val_Class %in% c("KPC_POSITIVE", "negative"), c("OTU")])
# which ones are shared by all groups?
test<-Reduce(intersect, list(all, val, cal))
	# combine these together (without duplicates)
otulist<-sort(unique(c(all, val, cal)))
	# now, get relevant information from lefse file
important_otus<-lf[lf$OTU %in% otulist, c("OTU", "taxname", "All_LDA", "All_pValue", "Cal_LDA", "Cal_pValue", "Val_LDA", "val_pValue")]
important_otus$genus<-gsub("Otu.*?_", "", important_otus$taxname)
row.names(important_otus)<-gsub("000", "", important_otus$OTU)
important_otus[3:8]<-signif(important_otus[3:8], 2)
important_otus[is.na(important_otus)]<-"ns"
important_otus<- important_otus[seq(dim(important_otus)[1],1),]

# step 2: create 2 dataframes for each kpc status group from shared file:
rawshared<-read.table(file="manuscript/Git/Data/rushCX_filtered.shared.frac.txt", sep="\t", header=TRUE, row.names=1)
cx<-read.table(file="manuscript/Git/Data/rushfinalCX_summary.w.pam.txt", header=TRUE, sep="\t")
fshared<-rawshared[, otulist ]
neglist<-as.character(cx[cx$cre_status %in% c("negative"), c("group")])
poslist<-as.character(cx[cx$cre_status %in% c("KPC_POSITIVE"), c("group")])
kpc_shared<-fshared[row.names(fshared) %in% poslist, ]
neg_shared<-fshared[row.names(fshared) %in% neglist, ]

# note: if you want logs for easier viewing:
log_kpc_shared <- log10(kpc_shared + 1)
log_neg_shared <- log10(neg_shared + 1)
log_kpc_shared<-log_kpc_shared[,order(colnames(log_kpc_shared),decreasing=TRUE)]
log_neg_shared<-log_neg_shared[,order(colnames(log_neg_shared),decreasing=TRUE)]

# step 3: plot
# OTU abundance differences: LOG transformed
par(mar=c(3,7,1,1), xaxs='r', mgp=c(2,1,0))
maxAb<-round(max(log_neg_shared), digits=2)
plot(1, type='n', ylim=c(0.8, (ncol(log_neg_shared)*2)-0.8), xlim=c(0,maxAb), 
     ylab='', xlab='log Relative Abundance', xaxt='n', yaxt='n', cex.lab=1)
index <- 1
for(i in colnames(log_neg_shared)){
  stripchart(at=index+0.35, log_neg_shared[,i], 
             pch=21, bg='midnightblue', method='jitter', jitter=0.15, cex=1, lwd=0.5, add=TRUE)
  stripchart(at=index-0.35, log_kpc_shared[,i], 
             pch=21, bg='gold', method='jitter', jitter=0.15, cex=1, lwd=0.5, add=TRUE)
  if (i != colnames(log_neg_shared)[length(colnames(log_neg_shared))]){
    abline(h=index+1, lty=2)
  }
  segments(mean(log_neg_shared[,i]), index+0.9, mean(log_neg_shared[,i]), index+0.1, lwd=2.5, col="black") #adds line for median
  segments(mean(log_kpc_shared[,i]), index-0.9, mean(log_kpc_shared[,i]), index-0.1, lwd=2.5, col="black")
  index <- index + 2
}
axis(side=1, at=c(0,maxAb), label=c(0,maxAb), cex.axis=1, tck=-0.02)
minors <- c(0.1,0.28,0.44,0.58,0.7,0.8,0.88,0.94,0.98)
axis(side=1, at=minors, label=rep('',length(minors)), tck=-0.01)
#axis(side=1, at=minors+1, label=rep('',length(minors)), tck=-0.01)
#axis(side=1, at=minors+2, label=rep('',length(minors)), tck=-0.01)
#axis(side=1, at=minors+3, label=rep('',length(minors)), tck=-0.01)
legend('bottomright', legend=c('negative', 'KPC(+)'),
       pch=c(21, 21), pt.bg=c('midnightblue','gold'), bg='white', pt.cex=1, cex=0.8)
axis(2, at=seq(1,index-2,2)+0.6, labels=rownames(important_otus), las=1, line=-0.5, tick=F, cex.axis=0.75)
formatted_taxa <- lapply(1:nrow(important_otus), function(x) bquote(paste(italic(.(important_otus$genus[x])), sep='')))
axis(2, at=seq(1,index-2,2), labels=do.call(expression, formatted_taxa), las=1, line=-0.5, tick=F, cex.axis=0.75, font=3) 
#italic_p <- lapply(1:length(important_otus$All_LDA), function(x) bquote(paste(italic('p = '), .(signif(important_otus$All_pValue[x], 2)), sep=' ')))
italic_p <- lapply(1:length(important_otus$All_LDA), function(x) bquote(paste(italic('p = '), .(important_otus$All_pValue[x]), sep=' ')))
axis(2, at=seq(1,index-2,2)-0.6, labels=do.call(expression, italic_p), las=1, line=-0.5, tick=F, cex.axis=0.5, font=3)
#mtext('B', side=2, line=2, las=2, adj=13, padj=-13, cex=1.7)

	
# step 4: let's add a 'heatmap' indicating the LDA value for each of the comparisons...
# this is a bit more complicated
important_otus<-lf[lf$OTU %in% otulist, c("OTU", "taxname", "All_LDA", "All_pValue", "Cal_LDA", "Cal_pValue", "Val_LDA", "val_pValue")]
important_otus$genus<-gsub("Otu.*?_", "", important_otus$taxname)
row.names(important_otus)<-gsub("000", "", important_otus$OTU)
otu.matrix<-as.matrix(important_otus[, c("All_LDA", "Cal_LDA", "Val_LDA")])

dat<-otu.matrix
my_palette <- colorRampPalette(c("yellow", "orange", "red")) (n=6)
breaks <- seq(min(dat, na.rm = T), max(dat, na.rm = T), length.out = 7)
heatmap.2(dat, trace="none", na.color = "grey97", scale="none", col = my_palette, breaks=breaks, 
	dendrogram="none", Colv=F, Rowv=F, density.info="none", 
	cexRow=0.8, cexCol=0.8, labCol=c("All", "Calibration", "Validation"))
	
## note: the heatmap call does not allow for a multi-figure panels, so just combined them in a vector-based editing program (Adobe illustrator)


```