### R script to plot bins
### arg1 = map_and_bin folder full path
### arg2 = metadata file full path

library(readr)
source("/work/apps/FastPCA.R")
library(Rtsne)
library(viridis)
library(fastcluster)
library(ape)
library(dplyr)
library(tibble)
library(vegan)
library(gplots)
library(readxl)
library(tidyr)
library(leaflet)
library(scales)

source("/data/dana/val2col.R")

args <- commandArgs(trailingOnly = TRUE)

prjdir = args[1]
setwd(file.path(prjdir))

metadata_file <- args[2]
meta <- read_excel(metadata_file, sheet = "nanopore log")

dir.create(file.path(prjdir,"plots"))
setwd(file.path(prjdir,"plots"))

pdf(file="plots.pdf", width=24, height=24)

#meta <- read_csv(metadata_file)
meta <- meta[!is.na(meta$flowcell) & !is.na(meta$barcode), ]
project = 'CMO2025'
#if (!(project == "."))  meta <- meta %>% filter(project == !!project)

meta$longitude = meta$longitude + sample(1:1000, nrow(meta)) / 1e6 #add jitter so we can zoom and see differences
meta$latitude = meta$latitude + sample(1:1000, nrow(meta)) / 1e6 #add jitter so we can zoom and see differences
meta$latitude[is.na(meta$latitude)] = 90
fakelong = -100 + as.numeric(as.factor(meta$flowcell)) + as.numeric(as.factor(meta$barcode)) /
  max(as.numeric(as.factor(meta$barcode))) / 3
meta$longitude[is.na(meta$longitude)] = fakelong[is.na(meta$longitude)]


depth.in = read_tsv(file.path(prjdir,"/bamdir/depths_jgi.txt"))
rownames(depth.in) = depth.in$contigName

mat.depth = depth.in[,1:3]
depth = depth.in[,seq(4,ncol(depth.in),2)]
depth.sm = depth[,colSums(depth>0)>1000]

depth.pca = FastPCA(depth^.25, 5) #50
rownames(depth.pca$x) = mat.depth$contigName
plot(depth.pca$x[,1:2],pch='.')

#depth.tsne = Rtsne(depth.pca$x, dims = 2, check_duplicates = F, pca = F, max_iter = 2500, verbose = T, num_threads = 8)
#rownames(depth.tsne$Y) = depth.in$contigName

#plot(depth.tsne$Y[,1:2],pch='.')

tnf.hd = read_tsv(file.path(prjdir,"tetra/tnfs.txt"))
tnfs = read_tsv(file.path(prjdir,"tetra/assembly.lrn"),col_names=F)
colnames(tnfs) = colnames(tnf.hd)
tnfs.mat = data.frame("contigName"=tnfs$seqid, "tnfs"=1)
rownames(tnfs) = tnfs.mat$contigName
tnfs$seqid = NULL

tnfs.pca = FastPCA(tnfs, 50)
plot(tnfs.pca$x[,1:2],pch='.')
plot(tnfs.pca$x[,1:2],pch='.',cex=rescale(rowSums(depth),c(0,3)))

rownames(tnfs.pca$x) = tnfs.mat$contigName

if(file.exists("tnfs.tsne.rds")) {
tnfs.tsne = readRDS("tnfs.tsne.rds")
} else {
tnfs.tsne = Rtsne(tnfs.pca$x, dims = 2, check_duplicates = F, pca = F, max_iter = 2500, verbose = T, num_threads = 8)
}
plot(tnfs.tsne$Y[,1:2],pch='.')
rownames(tnfs.tsne$Y) = tnfs.mat$contigName
saveRDS(tnfs.tsne,file="tnfs.tsne.rds")

bins.in = read_tsv(file.path(prjdir,"dastool_DASTool_bins/dastool_DASTool_contig2bin.tsv"),col_names=F)
colnames(bins.in) = c("contigName","dastoolbin")
bins.in$dastoolbin = paste0("B_",bins.in$dastoolbin)

bins.info = read_delim(file.path(prjdir,"dastool_DASTool_bins/kaiju.allbins.summary.tsv"),col_names=F,delim = "\t")
colnames(bins.info) = c("contigName","length","binid","classified","c2","taxid","taxstring")

#x=bins.info[rank(bins.info$length,ties="random"),]
bins.info.na = bins.info[!is.na(bins.info$taxstring),]
bins.info.na = bins.info.na[bins.info.na$taxstring != "Bacteria; NA; NA; NA; NA; NA; bacterium; ",]
x=bins.info.na[rev(rank(bins.info.na$length,ties="random")),]
y = x[!duplicated(x$binid),]
z = y$taxstring
names(z) = paste0("B_",y$binid)

mat = left_join(left_join(left_join(mat.depth,bins.info),tnfs.mat),bins.in)

mat$dastoolbin[is.na(mat$dastoolbin)] = 0
mat$binfactor = as.factor(mat$dastoolbin)
nmax = length(unique(mat$binfactor))
mat$cols = rainbow(nmax,alpha = 0.5)[as.numeric(mat$binfactor)]

cols.depth = mat$cols[!is.na(mat$dastoolbin)]
names.depth = mat$contigName[!is.na(mat$dastoolbin)]
length.depth = mat$contigLen[!is.na(mat$dastoolbin)]

# skipped depth tsne
#binned.depth = depth.tsne$Y[!is.na(mat$dastoolbin),]
#plot(depth.tsne$Y[,1:2],pch='.')
#points(binned.depth,cex=sqrt(length.depth)/100,pch=21,bg=cols.depth,col=NULL)

# mat$tnfrow = match(mat$contigName, tnfs.mat$contigName)
# pickrows = mat$tnfrow[(!is.na(mat$dastoolbin) & !is.na(mat$tnfs))]
# binned.tnfs = tnfs.tsne$Y[mat$contigName[pickrows],]
# cols.tnfs = mat[mat$contigName[pickrows],"cols"]
# names.tnfs = mat[mat$contigName[pickrows],"contigName"]
# length.tnfs = mat[mat$contigName[pickrows],"contigLen"]
# 
# plot(tnfs.tsne$Y[,1:2],pch='.')
# points(binned.tnfs,cex=sqrt(length.tnfs)/100,pch=21,bg=cols.tnfs,col=NULL)



bins_lut <- bins.in %>%
  mutate(bin = as.factor(dastoolbin)) %>%
  select(contigName, bin)

# build a color palette per bin
bin_levels <- levels(bins_lut$bin)
bin_cols <- setNames(rainbow(length(bin_levels), alpha = 0.8), bin_levels)

idx <- match(rownames(tnfs.tsne$Y), bins_lut$contigName)
bin_for_row <- bins_lut$bin[idx]            # factor with NAs for unbinned contigs
col_for_row <- bin_cols[as.character(bin_for_row)]
col_for_row[is.na(col_for_row)] <- "#00000020"  # faint grey for unbinned

size_for_row <- 1
if ("contigLen" %in% colnames(depth.in)) {
  len_idx <- match(rownames(tnfs.tsne$Y), depth.in$contigName)
  clen <- depth.in$contigLen[len_idx]
  # sane size scaling
  size_for_row <- pmax(0.4, sqrt(pmax(clen, 1))/300)
}
plot(tnfs.tsne$Y[,1], tnfs.tsne$Y[,2], 
     pch = 21, bg = col_for_row, col = NA,
     cex = size_for_row, 
     xlab = "t-SNE 1 (TNF)", ylab = "t-SNE 2 (TNF)",
     main = "TNF t-SNE colored by DAS Tool bin")

### hclust on samples
depth.min = depth[rowSums(depth > 0) > 4, colSums(depth > 0) > 4]
depth.sample.mat = t(prop.table(as.matrix(depth.min), margin = 2))
topk = min(50,nrow(depth.sample.mat))
depth.sample.pca = FastPCA(depth.sample.mat^.25,top.k = topk,center.use = TRUE)

depth.sample.dist = vegdist(depth.sample.mat,method="bray")
depth.sample.hclust = hclust(depth.sample.dist, method="ward.D")


fcbc <- strcapture(
  pattern = "^[^_]+_[^_]+_[^_]+_([^_]+)_[^_]+_(barcode\\d+)",
  x = depth.sample.hclust$labels,
  proto = list(flowcell = character(), barcode = character())
)
fcbc = left_join(fcbc, meta)
rownames(fcbc) = depth.sample.hclust$labels

rowlab = fcbc[rownames(depth.sample.mat), "sampleid"]
plot(depth.sample.hclust,labels = rowlab)

depth.sample.small = depth.sample.mat[,order(colMeans(depth.sample.mat),decreasing = TRUE)[1:3000]]
#depth.sample.small = depth.sample.mat[,max.col((depth.sample.mat))]
dim(depth.sample.small)

heatmap.2(depth.sample.small^.25,trace="none", col=viridis,
          # hclustfun=function(x) hclust(x, method="ward.D"),
          # distfun=function(x) vegdist(x, method="bray"),
          labRow = rowlab,
#          dendrogram = "row",
          Rowv = as.dendrogram(depth.sample.hclust)
          )

nclus.row = 9
sample.clusters = cutree(depth.sample.hclust,k=nclus.row)

heatmap.2(depth.sample.small^.25,trace="none", col=viridis,
          # hclustfun=function(x) hclust(x, method="ward.D"),
          # distfun=function(x) vegdist(x, method="bray"),
          labRow = rowlab,
          RowSideColors = rainbow(nclus.row)[as.numeric(sample.clusters)],
          #          dendrogram = "row",
          Rowv = as.dendrogram(depth.sample.hclust)
)



# aggregate by bin
depth.bins.mat = depth[match(bins.in$contigName,mat.depth$contigName),intersect(rownames(fcbc),colnames(depth))]
depth.bins.agg <- rowsum(as.matrix(depth.bins.mat), group = bins.in$dastoolbin)
depth.bins.agg = depth.bins.agg[rowSums(depth.bins.agg>0)>0,colSums(depth.bins.agg>0)>0]
dim(depth.bins.agg)
depth.bins = t(depth.bins.agg)
depth.bins.prop = prop.table(depth.bins.agg,margin=2)

depth.bins.dist = vegdist(depth.bins.prop^.25,method="bray")
depth.bins.hclust = hclust(depth.bins.dist, method="ward.D")
plot(depth.bins.hclust)

nclus.row=9
depth.bins.clusters = cutree(depth.bins.hclust,k=nclus.row)
rsc = rainbow(nclus.row)[as.numeric(depth.bins.clusters)]
names(rsc) = depth.bins.hclust$labels

pdf(file="hm.pdf", width=24, height=36)
collab = z[rownames(depth.bins.prop)]
rowlab = fcbc[colnames(depth.bins.prop), "sampleid"]
h1=heatmap.2(t(depth.bins.prop)^.5,trace="none", col=viridis,margins=c(100,20),
          hclustfun=function(x) hclust(x, method="ward.D"),
          distfun=function(x) vegdist(x, method="bray"),
          labRow = rowlab,
          labCol = collab,
          ColSideColors = rsc,
          Colv = as.dendrogram(depth.bins.hclust)
)

collab = z[rownames(depth.bins.prop)]
rowlab = fcbc[colnames(depth.bins.prop), "sampleid"]
h1=heatmap.2((depth.bins.prop)^.25,trace="none", col=viridis,margins=c(20,100),
          hclustfun=function(x) hclust(x, method="ward.D"),
          distfun=function(x) vegdist(x, method="bray"),
          labRow = collab, cexRow = 1,
          labCol = rowlab, cexCol = 1,
          RowSideColors = rsc,
          Rowv = as.dendrogram(depth.bins.hclust)
)
dev.off()


topk = min(25,ncol(depth.bins.prop))
depth.bins.pca = FastPCA(depth.bins.prop^.25,top.k = topk)
plot(depth.bins.pca$x,bg=rsc,pch=21)

fcbc$clus = depth.bins.clusters[rownames(fcbc)]
fcbc$color = rsc[rownames(fcbc)]
fcbc$size = round(depth.bins.prop[1,rownames(fcbc)]*100,1)

leaflet(fcbc) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(
    ~ longitude,
    ~ latitude,
    fillColor = ~color,
    color = ~color,
    fillOpacity = 0.3,
    radius = ~size,
    popup = paste(fcbc$location)
  ) 

#%>%
  # saveWidget(file = file.path(output_folder, file_name),
  #            selfcontained = FALSE)


# map for each bin sized to relabund
tnfs.in.depth = match(tnfs.mat$contigName,mat.depth$contigName)

dev.off()


# tsne for each sample sized to relabund
#pdf(file="samples.tsne.pdf",width=8,height=8)
png(file="samples.tsne%03d.png",width=1080,1080, units="px",bg="black")
op <- par(bg = "black",           # device/figure background
          mar=c(2.5, 0.5, 3, .5))  



for(samp in depth.sample.hclust$labels[depth.sample.hclust$order]) {
  if(! samp %in% colnames(depth.sm)) next
  print(samp)  
  plot(tnfs.tsne$Y[,1:2],pch='.',cex=0.1,col='gray15', 
       axes = FALSE, #ann = FALSE,
       xaxt = "n", yaxt = "n",
       xaxs = "i", yaxs = "i",      # no 4% padding
       bty  = "n",                  # no box
       main=fcbc[samp,"sampleid"], col.main="white"
  )

  size.abund = 0+(depth.sm[tnfs.in.depth, samp][[1]])^.25
  col.abund = val2col(size.abund,"viridis")
  points(tnfs.tsne$Y[,1:2], pch=21, bg=col.abund, col=NULL,cex=0.45*size.abund)

  
}
dev.off()



# tsne for each bin sized to relabund
#not working

#pdf(file="bins.tsne.pdf",width=8,height=8)
#for(bin in 1:ncol(depth.bins.prop)) {
#
#  pick.contigs = match(bins.in$contigName[bins.in$dastoolbin==colnames(depth.bins.prop)[bin]],tnfs.mat$contigName)
#  print(sum(pick.contigs))
#
#  plot(tnfs.tsne$Y[,1:2],pch='.')
#  points(tnfs.tsne$Y[pick.contigs,1:2],pch=21,col="blue",bg="blue",cex=2)
#
#}
#dev.off()


# idx <- match(rownames(depth.tsne$Y), bins_lut$contigName)
# bin_for_row <- bins_lut$bin[idx]            # factor with NAs for unbinned contigs
# col_for_row <- bin_cols[as.character(bin_for_row)]
# col_for_row[is.na(col_for_row)] <- "#00000020"  # faint grey for unbinned
# size_for_row <- 1
# if ("contigLen" %in% colnames(depth.in)) {
#   len_idx <- match(rownames(depth.tsne$Y), depth.in$contigName)
#   clen <- depth.in$contigLen[len_idx]
#   # sane size scaling
#   size_for_row <- pmax(0.4, sqrt(pmax(clen, 1))/1000)
# }

# plot(depth.tsne$Y[,1], depth.tsne$Y[,2], pch = 21, bg = col_for_row, col = NA,
#      cex = size_for_row, xlab = "t-SNE 1 (COV)", ylab = "t-SNE 2 (COV)",
#      main = "COV t-SNE colored by DAS Tool bin")

#comb.mat = merge(tnfs.mat, mat.depth, by=0, all.x=TRUE) %>% column_to_rownames("Row.names")

#comb.pca = merge(tnfs.pca$x, depth.pca$x, by=0, all.x=TRUE) %>% column_to_rownames("Row.names")
#comb.tsne = Rtsne(comb.mat, dims = 2, check_duplicates = F, pca = F, max_iter = 2500, verbose = T, num_threads = 8)
#rownames(comb.tsne$Y) = rownames(comb.pca)
#plot(comb.tsne$Y,pch='.')

# build a color palette per bin
# bin_levels <- levels(bins_lut$bin)
# bin_cols <- setNames(rainbow(length(bin_levels), alpha = 0.8), bin_levels)
# 
# idx <- match(rownames(comb.tsne$Y), bins_lut$contigName)
# bin_for_row <- bins_lut$bin[idx]            # factor with NAs for unbinned contigs
# col_for_row <- bin_cols[as.character(bin_for_row)]
# col_for_row[is.na(col_for_row)] <- "#00000020"  # faint grey for unbinned
# 
# size_for_row <- 1
# if ("contigLen" %in% colnames(depth.in)) {
#   len_idx <- match(rownames(comb.tsne$Y), depth.in$contigName)
#   clen <- depth.in$contigLen[len_idx]
#   # sane size scaling
#   size_for_row <- pmax(0.4, sqrt(pmax(clen, 1))/1000)
# }
# plot(comb.tsne$Y[,1], comb.tsne$Y[,2], pch = 21, bg = col_for_row, col = NA,
#      cex = size_for_row, xlab = "t-SNE 1 (TNF)", ylab = "t-SNE 2 (TNF)",
#      main = "TNF t-SNE colored by DAS Tool bin")


# comb.dist = dist(comb.tsne$Y)
# comb.tre = hclust(comb.tsne$Y, method = "single")


