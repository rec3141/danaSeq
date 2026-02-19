library(Rtsne)
source("~/apps/FastPCA.R")
library(readr)
library(scales)

args = commandArgs(trailingOnly=TRUE)

setwd(args[1])

tnfs = as.character(read.table("tnfs.txt",stringsAsFactors=F))
tnfs = tnfs[2:length(tnfs)]

lrn = read_table("Tetra_all.lrn",comment="%",col_names=F)
lrn$X1 = NULL

ATs = rowSums(lrn[,!grepl("G|C",tnfs)])
GCs = rowSums(lrn[,!grepl("A|T",tnfs)])
GCpct = GCs/(ATs+GCs)

lrn.pca = FastPCA(lrn,50)

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))

#This adds a column of color values
# based on the y values
GCcol <- rbPal(100)[as.numeric(cut(GCpct,breaks = 100))]
ptcol <- rbPal(100)[as.numeric(cut(lrn.pca$x[,1],breaks = 100))]
#ptcol <- rbPal(100)[as.numeric(cut(lrn.pca$x[,2],breaks = 100))]

#[sample(1:nrow(lrn),50000),]
#lrn.tsne = readRDS(file="tsne.RDS")
lrn.tsne = Rtsne(lrn.pca$x, check_duplicates=F, pca=F, verbose=T, num_threads=8)
saveRDS(file="tsne.RDS",lrn.tsne)

pdf(file="tsne.pdf", width=12, height=12)
plot(lrn.tsne$Y,pch='.', col=alpha(ptcol,0.1))
dev.off()

png(file="tsne.png", width=1200, height=1200)
plot(lrn.tsne$Y,pch=19, col=ptcol)
dev.off()

#tetids = read.table("Tetra_all.names",comment="%",row.names=1,stringsAsFactors=F)

