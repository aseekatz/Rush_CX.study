#### Rush CX study: Fig. 2, community diversity
##### Anna M. Seekatz
##### 12.18.17

###### Figures created:
	- Fig. 2A, B: diversity by KPC-Kp status, abx exposure
	- Fig. S1: diversity by different abx categories


####### Files used:
	- rushfinalCX_summary.w.pam.txt

######## Fig. 2, S1: Diversity

```
library(shape)

sums<-read.table(file="manuscript/Git/Data/rushfinalCX_summary.w.pam.txt", header=TRUE, sep="\t")

# graph diversity by CRE status:
group.col <- function(n) {
	colorvec <- vector(mode="character", length=length(n))
	for (i in 1:length(n)) {
	colorvec[i] = "light grey"
	if ( n[i] == "KPC_POSITIVE" ) {
	colorvec[i] = "gold"
	}
	if ( n[i] == "negative" ) {
	colorvec[i] = "midnightblue"
	}
	}
	c(colorvec)
	}
# Fig. 2A (by CRE status):
data<-droplevels(sums)
plot<-plot(invsimpson_03 ~ as.factor(cre_status), data = data, ylab=expression(paste("", lambda, " (diversity)")), 
	xlab="", xaxt="n", outline=FALSE, mpg=c(2,1,0), ylim=c(0,28), cex=1, cex.lab=1, cex.axis=1)
points(invsimpson_03~ jitter(as.numeric(cre_status, factor=0)), data = data, bg=group.col(data$cre_status), col="black", pch=21, cex=1)
names<-c("KPC(+)", "negative")
text(x =  seq(1,2,by=1), y = par("usr")[3]-1, srt = 45, adj = 1, labels = names, xpd = TRUE, cex=1)
#kruskal.test(invsimpson_03~cre_status, data=data)
wilcox.test(invsimpson_03~cre_status, data=data)
text(1.5,27, labels="Wilcoxon test, ns", cex=0.8)

### supplemental diversity figures

# let's look at the diversity by some different measures
colnames(data)
div.boxplot<- function (n, m) {
				plot(data[,n] ~ as.factor(data[,m]), data = data, ylab=names(data[n]), 
					xlab="", xaxt="n", outline=FALSE, mpg=c(2,1,0), ylim=c(0,max(data[,n])), cex=1, cex.lab=1, cex.axis=1)
				points(data[,n] ~ jitter(as.numeric(data[,m], factor=0)), data = data, bg=c("deepskyblue", "dodgerblue4"), col="black", pch=21, cex=1)
				print(as.character(unique(data[,m])))
				names<-as.character(unique(data[,m]))
				text(x =  seq(1,length(unique(data[,m])),by=1), y = par("usr")[3]-max(data[n], na.rm=TRUE)/11, srt = 45, adj = 1, labels = names, xpd = TRUE, cex=1)
				stat<-kruskal.test(data[,n]~data[,m], data)
				print(stat$p.value)
				text(length(unique(data[,m]))/1.3, max(data[,n]), labels=paste0("Kruskal-Wallis test \np < ", signif(stat$p.value, digits=3)), cex=0.5)
				mtext(names(data[m]), 1, line=2)
				}
#measures<-c("cre_status", "age", "sex", "Antibiotic_all", "Class_Carbapenem", "vancomycin_IV", "Gut_anaerobes", "High_biliary", "High_stool")
par(mfrow=c(2,2))
data[23:74] <- lapply(data[23:74], factor)
lapply(c("invsimpson_03"), FUN=div.boxplot, m=c("cre_status", "age", "sex", "Antibiotic_all", "Class_Carbapenem", "vancomycin_IV", "Gut_anaerobes", "High_biliary", "High_stool"))
	#note: this is not being called out correctly
lapply(c("invsimpson_03"), FUN=div.boxplot, m=c("Gut_anaerobes"))
lapply(c("invsimpson_03"), FUN=div.boxplot, m=c("High_biliary"))
lapply(c("invsimpson_03"), FUN=div.boxplot, m=c("High_stool"))
title('Diversity measures', outer=TRUE, line=-2)

# redid code to add color:
abx.col <- function(n) {
	colorvec <- vector(mode="character", length=length(n))
	for (i in 1:length(n)) {
	colorvec[i] = "light grey"
	if ( n[i] == "1" ) {
	colorvec[i] = "deepskyblue"
	}
	if ( n[i] == "0" ) {
	colorvec[i] = "dodgerblue4"
	}
	}
	c(colorvec)
	}
abx.div.boxplot<- function (n, m) {
				plot(data[,n] ~ as.factor(data[,m]), data = data, ylab=names(data[n]), 
					xlab="", xaxt="n", outline=FALSE, mpg=c(2,1,0), ylim=c(0,max(data[,n])), cex=1, cex.lab=1, cex.axis=1)
				points(data[,n] ~ jitter(as.numeric(data[,m], factor=0)), data = data, bg=abx.col(data[,m]), col="black", pch=21, cex=1)
				print(as.character(unique(data[,m])))
				names<-as.character(unique(data[,m]))
				text(x =  seq(1,length(unique(data[,m])),by=1), y = par("usr")[3]-max(data[n], na.rm=TRUE)/11, srt = 45, adj = 1, labels = names, xpd = TRUE, cex=1)
				stat<-kruskal.test(data[,n]~data[,m], data)
				print(stat$p.value)
				text(length(unique(data[,m]))/1.3, max(data[,n]), labels=paste0("Kruskal-Wallis test \np < ", signif(stat$p.value, digits=3)), cex=0.5)
				mtext(names(data[m]), 1, line=2)
				}
#measures<-c("cre_status", "age", "sex", "Antibiotic_all", "Class_Carbapenem", "vancomycin_IV", "Gut_anaerobes", "High_biliary", "High_stool")
par(mfrow=c(2,2))
data[23:74] <- lapply(data[23:74], factor)
lapply(c("invsimpson_03"), FUN=div.boxplot, m=c("Antibiotic_all"))
	#note: this is not being called out correctly
lapply(c("invsimpson_03"), FUN=div.boxplot, m=c("Gut_anaerobes"))
lapply(c("invsimpson_03"), FUN=div.boxplot, m=c("High_biliary"))
lapply(c("invsimpson_03"), FUN=div.boxplot, m=c("High_stool"))
title('Diversity measures', outer=TRUE, line=-2)

### Fig 2, A/B:
par(mfrow=c(1,2), xpd=TRUE, mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))
# A: cre status:
	# note: change div.boxplot to note 'group.col' and 'wilcox.test' instead of 'abx.col' and 'kruskal.test'
data<-droplevels(sums)
cre.div.boxplot<- function (n, m) {
				plot(data[,n] ~ as.factor(data[,m]), data = data, ylab=names(data[n]), 
					xlab="", xaxt="n", outline=FALSE, mpg=c(2,1,0), ylim=c(0,max(data[,n])), cex=1, cex.lab=1, cex.axis=1)
				points(data[,n] ~ jitter(as.numeric(data[,m], factor=0)), data = data, bg=group.col(data[,m]), col="black", pch=21, cex=1)
				print(as.character(unique(data[,m])))
				names<-as.character(unique(levels(data[,m])))
				text(x =  seq(1,length(unique(data[,m])),by=1), y = par("usr")[3]-max(data[n], na.rm=TRUE)/11, srt = 45, adj = 1, labels = names, xpd = TRUE, cex=1)
				stat<-wilcox.test(data[,n]~data[,m], data)
				print(stat$p.value)
				text(length(unique(data[,m]))/1.3, max(data[,n]), labels=paste0("Wilcoxon test \np < ", signif(stat$p.value, digits=3)), cex=0.5)
				mtext(names(data[m]), 1, line=2)
				}
lapply(c("shannon_03"), FUN=cre.div.boxplot, m=c("cre_status"))
# B: abx status:
data[23:74] <- lapply(data[23:74], factor)
abx.div.boxplot<- function (n, m) {
				plot(data[,n] ~ as.factor(data[,m]), data = data, ylab=names(data[n]), 
					xlab="", xaxt="n", outline=FALSE, mpg=c(2,1,0), ylim=c(0,max(data[,n])), cex=1, cex.lab=1, cex.axis=1)
				points(data[,n] ~ jitter(as.numeric(data[,m], factor=0)), data = data, bg=abx.col(data[,m]), col="black", pch=21, cex=1)
				print(as.character(unique(data[,m])))
				names<-as.character(unique(levels(data[,m])))
				text(x =  seq(1,length(unique(data[,m])),by=1), y = par("usr")[3]-max(data[n], na.rm=TRUE)/11, srt = 45, adj = 1, labels = names, xpd = TRUE, cex=1)
				stat<-wilcox.test(data[,n]~data[,m], data)
				print(stat$p.value)
				text(length(unique(data[,m]))/1.3, max(data[,n]), labels=paste0("Wilcoxon test \np < ", signif(stat$p.value, digits=3)), cex=0.5)
				mtext(names(data[m]), 1, line=2)
				}
lapply(c("shannon_03"), FUN=abx.div.boxplot, m=c("Antibiotic_all"))

### Fig 1S:
par(mfrow=c(1,3), xpd=TRUE, mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))
data[23:74] <- lapply(data[23:74], factor)
lapply(c("shannon_03"), FUN=abx.div.boxplot, m=c("Gut_anaerobes"))
lapply(c("shannon_03"), FUN=abx.div.boxplot, m=c("High_biliary"))
lapply(c("shannon_03"), FUN=abx.div.boxplot, m=c("High_stool"))

```



```
### Fig. 3A, by CRE status:
par(mfrow=c(1,2))
sums<-read.table(file="manuscript/Git/Data/rushfinalCX_summary.w.pam.txt", header=TRUE, sep="\t")
group.col <- function(n) {
  colorvec <- vector(mode="character", length=length(n))
  for (i in 1:length(n)) {
    colorvec[i] = "light grey"
    if ( n[i] == "KPC_POSITIVE" ) {
      colorvec[i] = "gold"
    }
    if ( n[i] == "negative" ) {
      colorvec[i] = "midnightblue"
    }
  }
  c(colorvec)
}
data<-droplevels(sums)
plot1<-plot(data$pcoa03_axis1, data$pcoa03_axis2, ylab="Axis 2 (6.8%)", xlab="Axis1 (8.3%)", data = data, pch=21, col="black", bg=group.col(data$cre_status), cex=1, cex.lab=1, cex.axis=1, ylim=c(-0.6, 0.6), xlim=c(-0.6, 0.6))
legend("bottomleft",legend=c("KPC-Kp (+)", "KPC-Kp (--)"), col=c("gold", "midnightblue"), cex=0.8, pch=19)

## note: need to add AMOVA value here?
# also: will add biplot here

### Fig. 3B: by PAM clustering:
# note: PAM clustering code is in code chunk below (complete that first)
pam.col <- function(n) {
  colorvec <- vector(mode="character", length=length(n))
  for (i in 1:length(n)) {
    colorvec[i] = "light grey"
    if ( n[i] == "1" ) {
      colorvec[i] = "turquoise"
    }
    if ( n[i] == "2" ) {
      colorvec[i] = "yellow"
    }
    if ( n[i] == "3" ) {
      colorvec[i] = "royalblue"
    }
  }
  c(colorvec)
}
plot2<-plot(data$pcoa03_axis1, data$pcoa03_axis2, ylab="Axis 2 (6.8%)", xlab="Axis1 (8.3%)", data = data, pch=21, col="black", bg=pam.col(data$phylopam3_cluster), cex=1, cex.lab=1, cex.axis=1, ylim=c(-0.6, 0.6), xlim=c(-0.6, 0.6))
legend("bottomleft",legend=c("PAM cluster 1", "PAM cluster 2", "PAM cluster 3"), col=c("turquoise", "yellow", "royalblue"), cex=0.8, pch=19)

# add biplot info:
bidata<-read.table(file="biplot/rushfinal.CXsamples.0.03.0.03.filter.spearman.corr.axes", header=TRUE)
tax<-read.table(file="manuscript/Git/Data/rushfinal.03_taxnames_allOTUs.txt", header=TRUE)
data<-merge(bidata, tax, by="OTU")
	# filter out by a cutoff for axes 1 and 2; nseqs
df<-data[data$length > 0.5 & data$p.value < 0.05 | data$length > 0.5 & data$p.value.1 < 0.05, ]
	# looks good to me...

# arrows from significant biplot data:
Arrows(0, 0, x1=-0.523011, y1=0.502029, lty=1, arr.type="triangle") #otu1
Arrows(0, 0, x1=0.442972, y1=-0.366383, lty=1, arr.type="triangle") #otu4
Arrows(0, 0, x1=0.427019, y1=0.528320, lty=1, arr.type="triangle") #otu5
Arrows(0, 0, x1=0.511414, y1=-0.055904, lty=1, arr.type="triangle") #otu30
Arrows(0, 0, x1=0.517663, y1=-0.022780, lty=1, arr.type="triangle") #otu32
Arrows(0, 0, x1=0.523796, y1=-0.093529, lty=1, arr.type="triangle") #otu74
Arrows(0, 0, x1= 0.584662, y1=-0.004479, lty=1, arr.type="triangle") #otu87
Arrows(0, 0, x1=0.461388, y1=0.122442, lty=1, arr.type="triangle") # otu90
Arrows(0, 0, x1=0.499862, y1=0.109177, lty=1, arr.type="triangle") #otu151
Arrows(0, 0, x1=0.069817, y1=-0.245070, lty=1, arr.type="triangle") #otu2
Arrows(0, 0, x1=-0.220818, y1=0.199955, lty=1, arr.type="triangle") #otu2

text(-0.523011, 0.502029, label="OTU001 \nEnterococcus", cex=.5, pos=1)
text(0.442972, -0.366383, label="OTU004 \nBacteroides", cex=.5, pos=1)
text(0.427019, 0.528320, label="OTU005 \nAkkermansia", cex=.5, pos=1)
text(0.511414, -0.055904, label="OTU030 \nClostridium XIVa", cex=.5, pos=1)
text(0.517663, -0.022780, label="OTU032 \nSubdoligranulum", cex=.5, pos=1)
text(0.523796, -0.093529, label="OTU074 \nFlavonifractor", cex=.5, pos=1)
text( 0.584662, -0.004479, label="OTU087 \nOscillibacter", cex=.5, pos=1)
text(0.461388, 0.122442, label="OTU090 \nClostridium XIVa", cex=.5, pos=1)
text(0.499862, 0.109177, label="OTU151 \nClostridium IV", cex=.5, pos=1)
text(0.069817, -0.245070, label="OTU002 \nEscherichia/Shigella", cex=.5, pos=1)
text(-0.220818, 0.199955, label="OTU002 \nEnterobacteriaceae", cex=.5, pos=1)

```



######## Fig. 3C and S2: PAM clustering relative abundances

```
#### pam clustering of samples:

library(cluster)
library(vegan)
library(labdsv)

# code for jensen-shannon distance:
all.pairs.jensen.shannon <- function(d) {
  kl.divergence <- function(pq, midpoint) {
    sum(ifelse(pq == 0.0, 0.0, pq * log(pq / midpoint)))
  }
 
  jensen.shannon <- function(row1, row2) {
    midpoint <- 0.5 * (row1 + row2)
    kl.divergence(row1, midpoint) + kl.divergence(row2, midpoint)
  }
 
  d <- as.matrix(d)
  num.rows <- nrow(d)
  result <- matrix(0, nrow=num.rows, ncol=num.rows, dimnames=list(rownames(d), rownames(d)))
  for (j in 1:(num.rows - 1)) {
    for (i in (j + 1):num.rows) {
      distance <- jensen.shannon(c(d[i, ]), c(d[j, ]))
      result[i, j] <- distance
      result[j, i] <- distance
    }
  }
  as.dist(result)
 }

# run your data through the code:
#sub.js.phylo <- all.pairs.jensen.shannon(sub.phylos.df)

# Because of the variability of OTUs across the data set, we will calculate PAM clusters based on phylotype, not OTU
# you could also do this with OTUs (messier)
	
### Now that a phylotype file has been created, you can use this to do pam clustering:
all.genera<-read.table(file="manuscript/Git/Data/rushfinalCX_all.genera.txt", sep="\t", header=TRUE)
cx<-read.table(file="manuscript/Git/Data/rushfinalCX_summary.txt", header=TRUE, sep="\t")

# run through script, and test with silhouette score for best clustering:
rownames(all.genera)<-all.genera$taxon
phylo.dd.frac<-as.matrix(t(all.genera[, 4:ncol(all.genera)]))
js.phylo <- all.pairs.jensen.shannon(phylo.dd.frac)

postscript('manuscript/Git/Figures/test.rushfinalCX.phylo_pamtest.ps',width=8,height=8,horizontal=TRUE)
for (nc in 2:10) {
plot(pam(js.phylo, k=nc))
}
dev.off()
	# anywhere from 2-5 clusters looks ok
	# for now, let's go with 5, which has the best silhouette score

# run with 5 clusters:
phylotype5 <- pam(js.phylo, k=5)										#the k= , will give you however many clusters you are defining
phylo5<-as.data.frame(phylotype5$silinfo$widths[,1])   					#this will list the samples, their assigned cluster, and the silhouette info for each sample
phylo5$SEQ_NAME<-rownames(phylo5)
colnames(phylo5)[1] <- c("phylopam5_cluster")
conv.summary<-merge(cx, phylo5, by.x=c("group"), by.y=c("SEQ_NAME"))

# run with 3 clusters:
phylotype3 <- pam(js.phylo, k=3)										#the k= , will give you however many clusters you are defining
phylo3<-as.data.frame(phylotype3$silinfo$widths[,1])   					#this will list the samples, their assigned cluster, and the silhouette info for each sample
phylo3$SEQ_NAME<-rownames(phylo3)
colnames(phylo3)[1] <- c("phylopam3_cluster")
conv.summary<-merge(conv.summary, phylo3, by.x=c("group"), by.y=c("SEQ_NAME"))

#write.table(conv.summary, file="manuscript/Git/Data/rushfinalCX_summary.w.pam.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

### now, create heatmap with PAM clustering:
# read in OTU and summary files:

library(RColorBrewer)
library(gplots)
library(vegan)
library(plyr)

rushsum<-read.table(file='manuscript/Git/Data/rushfinalCX_summary.w.pam.txt', header=TRUE, sep="\t")
shared<-read.table(file="manuscript/Git/Data/rushCX_filtered.shared.frac.txt", header=TRUE, row.names=1)

### Code for heatmap by sample (not included in manuscript):
otu.matrix<-as.matrix(shared)
	
# get relative abundance of file instead of raw counts:
otu.rel<-otu.matrix*100

# to filter out OTUs that are present in less than x%:
otu.rel.max<-apply(otu.rel,2,max)
otu.rel.filtered<-otu.rel[,otu.rel.max>0.020]
	#median, if you wanted to do that...
#otu.rel.med<-apply(otu.rel,2,median)
#otu.rel.filtered.med<-otu.rel[,otu.rel.med>0.0001]
			
# heatmap of top 40 (?) OTUs:
topotus<- otu.rel.filtered[, order(-colSums(otu.rel.filtered))]	
#top30<-topotus[, 1:30]
top40<-topotus[, 1:40]

# add row metadata:

rushsum<- rushsum[order(rushsum$group),]									#remember that the order of this must equal the order of the matrix read into the heatmap command
samples<-rownames(top40)										#this is the list of sample names in the heatmap matrix
#rushsum2<-rushsum[rushsum$seqID %in% samples , ]							#eliminate any samples from the original meta file to get the exact list (since not doing a merge)
top40<-top40[order(rownames(top40)), ]							#must be in same order
cbind(rownames(top40), as.character(rushsum$group))									#check that this is right

# giving some options to color by:
myCol5 <- colorRampPalette(brewer.pal(5,"BrBG"))(5)			# only for more than 3 colors
myCol3<-c("gold", "midnightblue")
PamCol<-c("turquoise", "yellow", "royalblue")
rushsum$color<-mapvalues(rushsum$cre_status, from = c("KPC_POSITIVE", "negative"), to = c("gold", "midnightblue"))
rushsum$phylopam5_cluster<-as.factor(rushsum$phylopam5_cluster)
rushsum$color2<-mapvalues(rushsum$phylopam5_cluster, from = c("1", "2", "3", "4", "5"), to = myCol5)
rushsum$color2<-mapvalues(rushsum$phylopam3_cluster, from = c("1", "2", "3"), to = PamCol)


# add column metadata (taxonomy info):
tax<-read.table(file="manuscript/Git/Data/rushCX.03_taxnames.txt", header=TRUE)

topotu.list<-colnames(top40)												#define the OTUs you want to keep (those in the top 98%)
filtered.tax<-tax[tax$OTU %in% topotu.list , ]							#filter out the taxonomy metadata file
filtered.tax<-droplevels(filtered.tax)										#sometimes R keeps levels that were discarded--this gets rid of them completely

filtered.tax$phylum<-factor(filtered.tax$phylum, levels=c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria", "unknown", "Verrucomicrobia"))
tax.col<- filtered.tax[order(filtered.tax$phylum, -filtered.tax$Size) , ]
tax.col$color<-mapvalues(tax.col$phylum, from = c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria", "unknown", "Verrucomicrobia"), to = c("firebrick1", "green4", "dodgerblue1","gold", "grey47", "hotpink"))
tax.col2<-t(tax.col)
phyla.col<-tax.col2[6, ]

col.order<-as.character(tax.col2[1, ])								#convert the order of the OTUs into a list
top40<-top40[,col.order]											#order your matrix by the ordered list
rbind(colnames(top40), tax.col2)									#check order of the matrix and list to ensure that colors will be correct (they should match!)

# by group:
top40<-top40[order(rushsum$cre_status), ]
rushsum<-rushsum[order(rushsum$cre_status), ]
# OR by phylotype:
top40<-top40[order(rushsum$phylopam3_cluster), ]
rushsum<-rushsum[order(rushsum$phylopam3_cluster), ]

# graph--subset samples:
my.col <- colorRampPalette(c("aliceblue","lightskyblue", "deepskyblue", "dodgerblue2", "blue4"))(n = 25)			#you can put in as many colors as you want here
my.breaks = c(seq(0,0.001,length=6),  								#just remember: your breaks must always = 1 more than the number of colors you are defining!
               seq(0.001,0.01,length=6),
               seq(0.01,0.1,length=6),
               seq(0.10,0.5,length=6),
               seq(0.50,1,length=6))
# does R 3.3.0 make these differently?!
#edited<-c(my.breaks[1:50], my.breaks[52:101], my.breaks[103:152], my.breaks[154:203], my.breaks[205:254])
edited.breaks<-my.breaks[!duplicated(my.breaks)]

### Fig. S2:
heatmap.2(top40,
                    notecol="black",     
                    density.info="none",
                    key.xlab="",
                    key.title="",  
                    trace="none",         
                    margins =c(2,4),    
                    col=my.col,        
                    breaks=edited.breaks,    
                    key=FALSE,  
                  	dendrogram="none",    
                    Colv=F,
                    Rowv=F,
                    srtCol=45,
                    cexRow= 0.7,
                    cexCol = 0.8,
                    keysize=1,
                    lwid = c(1,5),
                    lhei = c(2,5),
                    labCol=tax.col$taxname,
                    labRow="",
					RowSideColors = as.character(rushsum$color),
					ColSideColors = as.character(phyla.col)
					)
legend("topleft",legend=c("KPC-Kp (+)", "KPC-Kp (-)"), col=c("gold", "midnightblue"), cex=0.8, pch=19)
#legend("topleft",legend=rep(1:4, 1), col=myCol4, cex=0.8, pch=19)
legend("topright", c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria", "Verrucomicrobia", "unclassified" ), 
	col=c("firebrick1", "green4", "dodgerblue1","gold", "hotpink", "grey47"), pch=15, cex=0.8)
legend("left",legend=c("0-0.1", "0.1-1", "1-10", "10-50", "50-100"), col=c("aliceblue","lightskyblue", "deepskyblue", "dodgerblue2", "blue4"), cex=0.8, pch=19)
	

### since this is rather chunky, let's calculate the average relative abundance by pam cluster, then make a heatmap out of that:

### Fig 3C:
# combine data sets to reflect pam clusterings:
rushsum<-read.table(file='manuscript/Git/Data/rushfinalCX_summary.w.pam.txt', header=TRUE, sep="\t")
shared<-read.table(file="manuscript/Git/Data/rushCX_filtered.shared.frac.txt", header=TRUE, row.names=1)
pams<-rushsum[, c("group", "phylopam3_cluster")]
data<-merge(pams, shared, by.x="group", by.y='row.names')
rownames(data)<-data$group
data$group<-NULL

# get average relative abundance by PAM clustering:
df<-aggregate(data[, 2:303], list(data$phylopam3_cluster), mean)		# donot include rareOTUs (last column)
rownames(df)<-df$Group.1
otu.rel.filtered<-as.matrix(df[,2:ncol(df)]*100)

# limit to top30 OTUs, and define cluster colors:
topotus<- otu.rel.filtered[, order(-colSums(otu.rel.filtered))]	
top30<-topotus[, 1:30]
# giving some options to color by:
PamCol<-c("turquoise", "yellow", "royalblue")

# add column metadata (taxonomy info):
tax<-read.table(file="manuscript/Git/Data/rushCX.03_taxnames.txt", header=TRUE)

topotu.list<-colnames(top30)												#define the OTUs you want to keep (those in the top 98%)
filtered.tax<-tax[tax$OTU %in% topotu.list , ]							#filter out the taxonomy metadata file
filtered.tax<-droplevels(filtered.tax)										#sometimes R keeps levels that were discarded--this gets rid of them completely

filtered.tax$phylum<-factor(filtered.tax$phylum, levels=c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria", "unknown", "Verrucomicrobia"))
tax.col<- filtered.tax[order(filtered.tax$phylum, -filtered.tax$Size) , ]
tax.col$color<-mapvalues(tax.col$phylum, from = c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria", "unknown", "Verrucomicrobia"), to = c("firebrick1", "green4", "dodgerblue1","gold", "grey47", "hotpink"))
tax.col2<-t(tax.col)
phyla.col<-tax.col2[6, ]

col.order<-as.character(tax.col2[1, ])								#convert the order of the OTUs into a list
top30<-top30[,col.order]											#order your matrix by the ordered list
rbind(colnames(top30), tax.col2)

# now, graph this heatmap
#lowcol <- colorRampPalette(c("aliceblue", "cadetblue1", "cadetblue3"))(n = 3)
#medcol <- colorRampPalette(c("lightskyblue1", "lightskyblue"))(n = 2)
#highcol <- colorRampPalette(c("skyblue3", "skyblue3"))(n = 2)
#my.col <- c(lowcol, medcol, highcol)
my.col <- colorRampPalette(c("aliceblue", "lightblue1", "lightblue2", "skyblue3", "royalblue3"))(n = 10)
my.breaks = c(0, 0.1, 1, 2, 5, 7, 10, 15, 20, 50, 100)
#edited.breaks<-my.breaks[!duplicated(my.breaks)]
heatmap.2(top30,
                    notecol="black",     
                    density.info="none",
                    key.xlab="",
                    key.title="",  
                    trace="none",         
                    margins =c(2,4),    
                    col=my.col,        
                    breaks=my.breaks,    
                    key=FALSE,  
                  	dendrogram="none",    
                    Colv=F,
                    Rowv=T,
                    srtCol=45,
                    cexRow= 0.7,
                    cexCol = 0.8,
                    keysize=1,
                    lwid = c(1,5),
                    lhei = c(2,5),
                    labCol=tax.col$taxname,
                    labRow="",
					RowSideColors = as.character(PamCol),
					ColSideColors = as.character(phyla.col)
					)
legend("topright", c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria", "Verrucomicrobia", "unclassified" ), 
	col=c("firebrick1", "green4", "dodgerblue1","gold", "hotpink", "grey47"), pch=15, cex=0.8)
legend("topleft", paste("<", as.character(my.breaks[-1]), sep=" "), col=my.col, pch=15, cex=0.8)

```
