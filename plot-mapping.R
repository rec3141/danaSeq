library(Rtsne)
source("~/apps/FastPCA.R")
library(readr)
library(scales)

mat.in=read.delim("assembly.fasta.depth.txt",row.names=1)
mat = mat.in[,seq(3,ncol(mat.in),2)]

tsne = Rtsne(mat, check_duplicates=F)

plot(tsne$Y[,1:2],col=apply(mat,1,which.max))



###

system("head -n4 Tetra_all_clean_1000.lrn | tail -n1 | cut -f2- -d' ' > tnfs.txt")
tnfs = as.character(read.table("tnfs.txt",stringsAsFactors=F))
tnfs = tnfs[2:length(tnfs)]

lrn = read_table("Tetra_all_clean_1000.lrn",comment="%",col_names=F)
lrn$X1 = NULL

ATs = rowSums(lrn[,!grepl("G|C",tnfs)])
GCs = rowSums(lrn[,!grepl("A|T",tnfs)])
GCpct = GCs/(ATs+GCs)

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))

#This adds a column of color values
# based on the y values
GCcol <- rbPal(100)[as.numeric(cut(GCpct,breaks = 100))]
ptcol <- rbPal(100)[as.numeric(cut(lrn.pca$x[,1],breaks = 100))]
ptcol <- rbPal(100)[as.numeric(cut(lrn.pca$x[,2],breaks = 100))]


lrn.pca = FastPCA(lrn,50)

lrn.tsne = Rtsne(lrn.pca$x, check_duplicates=F, pca=F, verbose=T, num_threads=8)

tetids = read.table("Tetra_all_clean_1000.names",comment="%",row.names=1,stringsAsFactors=F)
binmem = read.table("bin.MemberMatrix.txt",row.names=1,head=T) #which bin is a read in
binids = paste0("bin",binmem$ClusterId[match(tetids$V3,rownames(binmem))]) #which bin is a shredded read in
bincol = rainbow(length(unique(binids)))[as.numeric(as.factor(binids))] # color of shredded read bin

binmem.5k = read.table("bin5k.MemberMatrix.txt",row.names=1,head=T) #which bin is a read in
binids.5k = paste0("bin5k",binmem.5k$ClusterId[match(tetids$V3,rownames(binmem.5k))]) #which bin is a shredded read in
bincol.5k = rainbow(length(unique(binids.5k)))[as.numeric(as.factor(binids.5k))] # color of shredded read bin

depths.in = read_table("depths.txt")
depths = data.frame(depths.in[match(tetids$V3, depths.in$contigName),],stringsAsFactors=F)

depthcount = rowSums(depths[,seq(4,ncol(depths),2)]>0)
depthcount[depthcount<1] = NA
freqcol <- rbPal(max(depthcount,na.rm=T))[as.numeric(cut(depthcount,breaks = max(depthcount,na.rm=T)))]
plot(lrn.tsne$Y[,1:2],pch=19,col=freqcol, cex=0.1)

plot(lrn.tsne$Y[,1:2],pch=19,col=rgb(0,0,1,0.1), cex=0.1)
plot(lrn.tsne$Y[,1:2],pch=19,col=GCcol, cex=0.1)
plot(lrn.tsne$Y[,1:2],pch=19,col=ptcol, cex=0.1)



pdf(file="tsne.bins5k.pdf",width=8, height=8)
# plot(lrn.tsne$Y[,1:2],pch=21,bg=bincol,col=NULL,cex=0.1)
 plot(lrn.tsne$Y[,1:2],pch=21,bg=bincol.5k,col=NULL,cex=0.1)
#for(i in unique(binids)) {
for(i in unique(binids.5k)) {
print(i)
 plot(lrn.tsne$Y[,1:2], pch='.', col='black', cex=0.1,main=i)
# points(lrn.tsne$Y[binids %in% i,1:2],pch=21,bg=bincol[binids %in% i],col=NULL,cex=0.5)
 points(lrn.tsne$Y[binids.5k %in% i,1:2],pch=21,bg=bincol.5k[binids.5k %in% i],col=NULL,cex=0.2)
}
dev.off()

pdf(file="tsne.pdf",width=8, height=8)
ptcol <- rbPal(100)[as.numeric(cut(lrn.pca$x[,1],breaks = 100))]
plot(lrn.tsne$Y[,1:2],pch=19,col=alpha(ptcol,0.1), cex=0.1)
for(i in seq(4,ncol(depths),2)) {
plot(lrn.tsne$Y[,1:2],pch=19,col=rgb(0,0,1,0.1), cex=as.numeric(depths[,i])/4)
}
dev.off()


### color with metabat bins
#lrn = read_table("Tetra_assembly_1000.lrn",sep="\t",comment="%",row.names=1)
#lrn = read_table("Tetra_assembly_1000.lrn",comment="%")

#binmem = read.table("bin.MemberMatrix.txt",row.names=1,head=T)
#binids = paste0("bin",binmem$ClusterId[match(tetids$V3,rownames(binmem))])
#bincol = rainbow(length(unique(binids)))[as.numeric(as.factor(binids))]
#plot(lrn.tsne$Y[,1:2],pch=21,bg=bincol,col=NULL,cex=0.5)


sketchtax = read.table("./../sketch_bin.tax",sep="\t",stringsAsFactors=F)
rownames(sketchtax) = sketchtax$V1

sketchprottax = read.table("./../sketch_prot_bin.tax",sep="\t",stringsAsFactors=F)
rownames(sketchprottax) = sketchprottax$V1

bindepths.in = read_table("bin_depths.txt")
bindepths = data.frame(bindepths.in[match(tetids$V3, bindepths.in$contigName),seq(4,ncol(bindepths.in),2)],stringsAsFactors=F)

sketchcol = sapply(sketchtax$V1,function(x) grep(x,colnames(bindepths)))
sketchprotcol = sapply(sketchprottax$V1,function(x) grep(x,colnames(bindepths)))


## plot each bin in blue
pdf(file="tsne_bins.pdf",width=8, height=8)
ptcol <- rbPal(100)[as.numeric(cut(lrn.pca$x[,1],breaks = 100))]
plot(lrn.tsne$Y[,1:2],pch=19,col=alpha(ptcol,0.1), cex=0.1)

plot(lrn.tsne$Y[,1:2],pch=19,col=alpha(ptcol,0.1), cex=0.1)
points(lrn.tsne$Y[,1:2],pch=19,col=rgb(0,0,1,1), cex=as.numeric(rowSums(bindepths)))

for(i in order(colSums(bindepths),decreasing=T)) {
print(i)
tax = sketchtax[names(which(sketchcol==i)),"V2"]
prottax = sketchprottax[names(which(sketchprotcol==i)),"V2"]

plot(lrn.tsne$Y[,1:2],pch='.',col='black', cex=0.1, main=paste0(colnames(bindepths)[i],"\n",tax,"\n",prottax))
points(lrn.tsne$Y[,1:2],pch=19,col=rgb(0,0,1,1), cex=as.numeric(bindepths[,i]))
}
dev.off()


### plot full color map for each bin separately
ptcol <- rbPal(100)[as.numeric(cut(lrn.pca$x[,1],breaks = 100))]
for(i in 1:ncol(bindepths)) {
pdf(file=paste0("tsne_bins_",i,".pdf"),width=8, height=8)

print(i)
tax = sketchtax[names(which(sketchcol==i)),"V2"]
prottax = sketchprottax[names(which(sketchprotcol==i)),"V2"]

plot(lrn.tsne$Y[,1:2],pch=19,col=alpha(ptcol,0.1), cex=0.1,main=paste0(colnames(bindepths)[i],"\n",tax,"\n",prottax))
#plot(lrn.tsne$Y[,1:2],pch='.',col='black', cex=0.1, main=paste0(colnames(bindepths)[i],"\n",tax,"\n",prottax))
points(lrn.tsne$Y[,1:2],pch=19,col=rgb(0,1,0,1), cex=as.numeric(bindepths[,i]))
dev.off()
}


### color by PCA

for(i in 1:20) {
pdf(file=paste0("tsne_pca",i,".pdf"),width=8, height=8)
ptcol <- rbPal(100)[as.numeric(cut(lrn.pca$x[,i],breaks = 100))]
plot(lrn.tsne$Y[,1:2],pch=19,col=alpha(ptcol,0.1), cex=0.1)
dev.off()
}




### just treat the cleaned-up fasta file as a co-assembly and then map each read to that?
### then do TNF binning on that? is that much worse than assembling first for a quick view of diversity?

library(vegan)
plot(hclust(vegdist(t(depth.mat),method="bray"),method="ward"))


lrnfiles = paste0("./../",list.files(pattern="Tetra_bin.*.lrn",path="./../"))

pdf(file="binlrn.pdf")
for(file in lrnfiles) {
binlrn = read_table(file,comment="%",col_names=F)
binlrn$X1 = NULL
hist(apply(binlrn,2,sd),breaks=seq(0,0.025,0.001),main=file)
}
dev.off()




#### plotting samples
samdepths.in = read_table("./../samples_bins/depths.txt")
#samdepths = data.frame(samdepths.in[match(samtetids$V3, samdepths.in$contigName),],stringsAsFactors=F)
samdepths = samdepths.in[,seq(4,ncol(samdepths.in),2)]

system("grep '>' ./../samples_bins/bins.all.fa | cut -f2 -d'>' > ./../samples_bins/oldbinnames.txt")
obn = read.table(file="./../samples_bins/oldbinnames.txt",stringsAsFactors=F)
obn$contig = paste0("contig_",rownames(obn))
obn$bin = matrix(unlist(strsplit(obn$V1,split="_",fixed=T)),byrow=T,ncol=2)[,1]

samagg = aggregate(samdepths,by=list(obn$bin),mean)
rownames(samagg) = samagg$Group.1
samagg$Group.1 = NULL

heatmap(as.matrix(t(samagg))^.5,scale="row")




### plotting BLAST hits

blast.in = read.table("~/Desktop/blasting/myc_cluster_top.txt",stringsAsFactors=F)
pdf(file=paste0("tsne_blast.pdf"),width=8, height=8)
ptcol <- rbPal(100)[as.numeric(cut(lrn.pca$x[,1],breaks = 100))]
plot(lrn.tsne$Y[,1:2],pch=19,col=alpha(ptcol,0.1), cex=0.1)
points(lrn.tsne$Y[match(blast.in$V2,tetids$V2),1:2], pch=as.numeric(as.factor(blast.in$V1)))
dev.off()


### plotting BLAST hits per sample

blast.in = read.table("~/Desktop/blasting/myc_cluster_top.txt",stringsAsFactors=F)

pdf(file=paste0("tsne_blast.pdf"),width=8, height=8)
ptcol <- rbPal(100)[as.numeric(cut(lrn.pca$x[,1],breaks = 100))]
plot(lrn.tsne$Y[,1:2],pch=19,col=alpha(ptcol,0.1), cex=0.1)
points(lrn.tsne$Y[match(blast.in$V2,tetids$V2),1:2], pch=as.numeric(as.factor(blast.in$V1)))
dev.off()


 ## plot each sample
 # what happened to Cyano bin genes? they went away...
pdf(file="tsne_blast_sample.pdf",width=8, height=8)

depths.sub = depths[,seq(4,ncol(depths),2)]
ptcol <- rbPal(100)[as.numeric(cut(lrn.pca$x[,1],breaks = 100))]
pt.gene  = match(blast.in$V2, tetids$V2)
pt.pch = as.numeric(as.factor(blast.in$V1))

plot(lrn.tsne$Y[,1:2],pch=19,col=alpha(ptcol,0.1), cex=rowSums(depths.sub)/4)
pt.pick = which(rowSums(depths.sub[pt.gene,])>0)
points(lrn.tsne$Y[pt.gene[pt.pick],1:2], pch=pt.pch[pt.pick])

for(i in 1:ncol(depths.sub)) {
plot(lrn.tsne$Y[,1:2],pch=19,col=alpha(ptcol,0.1), cex=as.numeric(depths.sub[,i])/4, main=colnames(depths.sub)[i])
pt.pick = which(depths.sub[pt.gene,i]>0)

points(lrn.tsne$Y[pt.gene[pt.pick],1:2], pch=pt.pch[pt.pick])

}
dev.off()

## histogram of frequencies
#most have 1, some have up to 6 on one read

blast.ord = blast.in[order(blast.in$V2,blast.in$V1),]
blast.rle = unclass(rle(blast.ord$V2))

#> hist(rle(blast.ord$V2)$lengths,breaks=c(-1:6))$counts
#[1]   0 260  67 121  96   6  12

# probably whole pathway
pick.full = blast.rle$values[blast.rle$lengths>4]
plot(lrn.tsne$Y[,1:2],pch=19,col=alpha(ptcol,0.1), cex=rowSums(depths.sub)/4)
pt.pick = intersect(pt.gene, match(pick.full,tetids$V2))
points(lrn.tsne$Y[pt.pick,1:2], pch=21, col='black', bg=rgb(1,1,1,0.3), cex=2, lwd=2)

write.table(pick.full,file="myc_genes_full.txt",col.names=F,quote=F,row.names=F)
system("~/apps/bbmap/filterbyname.sh in=~/Desktop/out_data/all_clean.fasta include=t substring=t names=myc_genes_full.txt out=myc_genes_full.fasta")




### quadratic sampling for similarity estimation
xy.out = NULL

for(split in 1:5) {
xymax = ceiling(1.1*max(abs(range(lrn.tsne$Y))))

xcut = cut(lrn.tsne$Y[,1],breaks=seq(-1*xymax,xymax,length.out=1+2^split))
ycut = cut(lrn.tsne$Y[,2],breaks=seq(-1*xymax,xymax,length.out=1+2^split))
xycut = as.numeric(as.factor(paste0(xcut,"_",ycut)))

xyagg = aggregate(depths.sub,by=list(xycut),FUN=mean)

xy.out = cbind(xy.out,t(xyagg))

}

xy.out = xy.out[2:nrow(xy.out),]
xy.pca = FastPCA(xy.out,Df=50)




## plot each sample
#guanitoxin 
# the genes are not present on the same contig, they are not part of the BGC
# not unexpected, they searched thousands of datasets and only found it in a few
pdf(file="tsne_sample_guanitoxin.pdf",width=8, height=8)
#guan.in = read.table("~/Desktop/blasting/guanitoxin/guanitoxin.best.txt",stringsAsFactors=F)
#colnames(guan.in) = c("qseqid","sseqid","evalue","qcovs")

gnt.all = NULL
for(file in list.files(".","tblastn_Gnt.*best2.out")) {
gnt.in = read.table(file,stringsAsFactors=F)
gnt.all = rbind(gnt.all,cbind(gnt.in, substr(file,start=9,stop=12)))
}

colnames(gnt.all) = c("qseqid","sseqid","evalue","qcovs", "gene")
guan.in = gnt.all

#guan.ord = guan.in[order(guan.in$sseqid,guan.in$qseqid),]
guan.ord = guan.in[order(guan.in$sseqid,guan.in$gene),]
guan.rle = unclass(rle(guan.ord$sseqid))

#> hist(rle(guan.ord$V2)$lengths,breaks=c(-1:6))$counts

# probably whole pathway
guan.full = guan.rle$values[guan.rle$lengths>1]

depths.sub = depths[,seq(4,ncol(depths),2)]
ptcol <- rbPal(100)[as.numeric(cut(lrn.pca$x[,1],breaks = 100))]
#pt.gene  = match(guan.in$sseqid, tetids$V3) #match all
pt.gene  = match(guan.full, tetids$V3) #match best of best

#pt.pch = as.numeric(as.factor(guan.in$qseqid))
pt.pch = rep(21,nrow(guan.in))

plot(lrn.tsne$Y[,1:2],pch=19,col=alpha(ptcol,0.1), cex=rowSums(depths.sub)/4)
pt.pick = which(rowSums(depths.sub[pt.gene,])>0)
points(lrn.tsne$Y[pt.gene[pt.pick],1:2], pch=pt.pch[pt.pick])

for(i in 1:ncol(depths.sub)) {
plot(lrn.tsne$Y[,1:2],pch=19,col=alpha(ptcol,0.1), cex=as.numeric(depths.sub[,i])/4, main=colnames(depths.sub)[i])
pt.pick = which(depths.sub[pt.gene,i]>0)

points(lrn.tsne$Y[pt.gene[pt.pick],1:2], pch=pt.pch[pt.pick])

}
dev.off()

#


