library(Rtsne)
source("~/apps/FastPCA.R")
library(readr)
library(scales)
library(fastcluster)
# args = commandArgs(trailingOnly=TRUE)

# setwd(args[1])
setwd("/data/work/dana/out_data_project_CMO/")

tnfs = as.character(read.table("tnfs.txt",stringsAsFactors=F))
tnfs = tnfs[2:length(tnfs)]

lrn = read_table("tetra/Tetra_all_clean_1500_1500.lrn",comment="%",col_names=F)
lrn$X1 = NULL

ATs = rowSums(lrn[,!grepl("G|C",tnfs)])
GCs = rowSums(lrn[,!grepl("A|T",tnfs)])
GCpct = GCs/(ATs+GCs)

lrn.pca = FastPCA(lrn,50)

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))

#This adds a column of color values
# based on the y values
GCcol <- rbPal(100)[as.numeric(cut(GCpct,breaks = 100))] #GC coloring
ptcol <- rbPal(100)[as.numeric(cut(lrn.pca$x[,1],breaks = 100))] #PC1 coloring
#ptcol <- rbPal(100)[as.numeric(cut(lrn.pca$x[,2],breaks = 100))]

lrn.pca.sub = lrn.pca$x[sample(1:nrow(lrn),100000),]
lrn.tsne = Rtsne(lrn.pca.sub, check_duplicates=F, pca=F, verbose=T, num_threads=8, max_iter = 2500)
saveRDS(file="tsne.RDS",lrn.tsne)
lrn.tsne = readRDS(file="tsne.RDS")
rownames(lrn.tsne$x) = rownames

pdf(file="tsne.pdf", width=8, height=8)
plot(lrn.tsne$Y,pch='.', col=alpha(GCcol,0.1))
dev.off()

png(file="tsne.png", width=1200, height=1200)
plot(lrn.tsne$Y,pch=19, col=ptcol)
dev.off()



#### plot bins on tsne
tetids = read.table("tetra/Tetra_all_clean_1500_1500.names",comment="%",row.names=1,stringsAsFactors=F)
#binmem = read.table("metabat/bin.MemberMatrix.txt",row.names=1,head=T) #which bin is a read in
binmem = read.table("metabat/bin.contigs.txt",head=F) #which bin is a read in
binids = binmem$V1[match(tetids$V3,binmem$V2)] #which bin is a shredded read in
bincol = rainbow(length(unique(binids)))[as.numeric(as.factor(binids))] # color of shredded read bin

sketchtax = read.table("sketch/refseq_bin.tax",sep="\t",stringsAsFactors=F)
rownames(sketchtax) = sketchtax$V1

sketchprottax = read.table("sketch/prot_bin.tax",sep="\t",stringsAsFactors=F)
rownames(sketchprottax) = make.names(sketchprottax$V1,unique = TRUE)

sketchcol = sapply(sketchtax$V1,function(x) grep(x,binmem$V1))
sketchprotcol = sapply(sketchprottax$V1,function(x) grep(x,binmem$V1))

contiglist = lapply(sketchcol, FUN = function(x) binmem$V2[x])

# plot each bin in blue
pdf(file="tsne_bins.pdf",width=8, height=8)
ptcol <- rbPal(100)[as.numeric(cut(lrn.pca$x[,1],breaks = 100))]
plot(lrn.tsne$Y[,1:2],pch=19,col=alpha(ptcol,0.1), cex=0.1)
#plot(lrn.pca$x[,1:2],pch=19,col=alpha(ptcol,0.1), cex=0.1)

plot(lrn.tsne$Y[,1:2],pch=19,col=alpha(ptcol,0.1), cex=0.1)
# plot(lrn.pca$x[,1:2],pch=19,col=alpha(ptcol,0.1), cex=0.1)
#points(lrn.tsne$Y[,1:2],pch=19,col=rgb(0,0,1,1), cex=as.numeric(rowSums(bindepths)))
points(lrn.tsne$Y[,1:2],pch=19,col=rgb(0,0,1,1), cex=0.1)
# points(lrn.pca$x[,1:2],pch=19,col=rgb(0,0,1,1), cex=0.1)

for(i in names(contiglist)) {
  print(i)
  tax = sketchtax[(which(names(sketchcol)==i)),"V2"]
  # prottax = sketchprottax[(which(names(sketchprotcol)==i))[1],"V2"]
  prottax=sketchprottax[(which(names(sketchcol)==i)),"V2"]
  pick.point = rep(0,nrow(lrn))
  pick.point[match(contiglist[[i]],tetids$V3)] = 0.2
  plot(lrn.tsne$Y[,1:2],pch='.',col='grey', cex=0.1, main=paste0(i,"\n",tax,"\n",prottax))
  points(lrn.tsne$Y[,1:2],pch=19,col=ptcol, cex=pick.point)
  # plot(lrn.pca$x[,1:2],pch='.',col='black', cex=0.1, main=paste0(tax,"\n",prottax))
  # points(lrn.pca$x[,1:2],pch=19,col=rgb(0,0,1,1), cex=c(0,1)[pick.point])
}
dev.off()





