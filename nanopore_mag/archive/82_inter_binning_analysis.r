#install.packages(c("gatepoints","Rtsne","Matrix","rgl","scales"))
#install_github("zarquon42b/Morpho", local=FALSE)
# had to module load some things to get these to compile
# library(umap)
# library(NbClust)
# library(mclust)
# library(cluster)
# library(MCL)
# library("FactoMineR")
# library("factoextra")
# library(klaR)
# library(cba)
# library(parallelDist)

library(gatepoints)
library(Rtsne)
library(Morpho)
library(Matrix)
library(data.table)
library(gplots)
library(rgl)
library(scales)
library(parallel)
library(igraph)
source("bin/FastPCA.R")

prefix="scaffolds_gt2000"

foamdir="annotation/foam"
kaijudir="annotation/kaiju"

kh <- read.table(file.path(foamdir,"kegg.headers.txt"),head=T)
fh <- read.table(file.path(foamdir,"foam.headers.txt"),head=T)

foamfile <- paste0(prefix,".foam.depths.tsv")
keggfile <- paste0(prefix,".kegg.depths.tsv")

foamDT <- fread(file.path(kaijudir,foamfile),sep="\t",head=F,stringsAsFactors=F)
keggDT <- fread(file.path(kaijudir,keggfile),sep="\t",head=F,stringsAsFactors=F)

colnames(foamDT) <- colnames(fh)
colnames(keggDT) <- colnames(kh)

oldrows.foamDT <- make.unique(foamDT$foam_hmm)
oldrows.keggDT <- make.unique(keggDT$foam_hmm)

rownames(foamDT) <- make.unique(foamDT$contig)
rownames(keggDT) <- make.unique(keggDT$contig)

foam.meta <- as.data.frame(foamDT[,1:11,with=F])
kegg.meta <- as.data.frame(keggDT[,1:13,with=F])

rownames(foam.meta) <- rownames(foamDT)
rownames(kegg.meta) <- rownames(keggDT)

foam.num <- as.data.frame(foamDT[,12:ncol(foamDT),with=F])
foam.num <- foam.num[,grep(".var",colnames(foam.num),invert=T)]
rownames(foam.num) <- make.unique(foamDT$orf)

kegg.num <- as.data.frame(keggDT[,14:ncol(keggDT),with=F])
kegg.num <- kegg.num[,grep(".var",colnames(kegg.num),invert=T)]
# 
# rownames(kegg.num) <- rownames(keggDT)


### dimensional reduction
#first with Fast PCA
covtab <- unique(foam.num)
cov.pca <- FastPCA(covtab,50)
cov.pca.uniq <- cov.pca$x
rownames(cov.pca.uniq) <- rownames(covtab)

# then Barnes-Hut t-SNE in 3D
cov.tsne3d.list <- Rtsne(cov.pca.uniq, dims=3, max_iter=1000, verbose=TRUE, theta=0.5, pca=FALSE)
saveRDS(cov.tsne3d.list,file="annotation/bins/tsne-embed.rds")
# cov.tsne3d.list <- readRDS(file="annotation/bins/tsne-embed.rds")
cov.tsne3d <- cov.tsne3d.list$Y
rownames(cov.tsne3d) <- rownames(covtab)

# stop()

# save.mat <- Matrix(0,nrow=nrow(covtab),ncol=nrow(covtab),sparse=T)
# rownames(save.mat) <- rownames(covtab)
# colnames(save.mat) <- rownames(covtab)

# save.class <- matrix(0,nrow=nrow(covtab),ncol=1000)
# rownames(save.class) <- rownames(covtab)

seqs <- rev(seq(40,120,1))
seqs <- 40:50 #have to have different numbers or hack below to find clusters will fail
save.class <- matrix(0,nrow=nrow(covtab),ncol=max(seqs))
rownames(save.class) <- sort(rownames(covtab))
for (count in seqs) {
	print(count)
	dms <- count

	fkm.tsne3d <- fastKmeans(cov.tsne3d, dms, iter.max = 10000, project = TRUE, threads = 0)
	plot3d(cov.tsne3d,pch='.',bg=NULL,col=alpha("grey",0.4),size=0.2,main=count,xlim=range(cov.tsne3d[,1]),ylim=range(cov.tsne3d[,2]),zlim=range(cov.tsne3d[,3]),axes=F)
	plot3d(cov.tsne3d,pch=21,col=alpha(fkm.tsne3d$class,0.5),cex=4,main=count,add=T)
# 	for(m in 1:dms) {
	save.tmp <- mclapply(X=1:dms,mc.cores=14,mc.preschedule=F, 
		FUN=function(m) {
		picked <- (fkm.tsne3d$class==m)
# 		hc <- hclust(dist(cov.tsne3d[picked,]))
 		hc <- hclust(dist(cov.tsne3d[picked,]),method="ward.D2")
		ct <- cutree(hc,h=max(hc$height)/5)
# 		plot(hc,labels=F)
# 	plot3d(cov.tsne3d[picked,],pch=21,col=alpha(ct,0.5),cex=4,main=count,add=T)
		})
		
	save.cls <- lapply(1:dms,function(m) paste(count,m,save.tmp[[m]],sep="_"))
	save.cls <- do.call(c,save.cls)
	name.cls <- lapply(1:dms,function(m) names(save.tmp[[m]]))
	name.cls <- do.call(c,name.cls)
	names(save.cls) <- name.cls
	save.class[,count] <- save.cls[order(names(save.cls))]
}

save.class <- save.class[,seqs]
save.test <- save.class[match(rownames(covtab),rownames(save.class)),]

# for some reason gives cluster1 too many...
# kmodes.out <- mclapply(X=seq(50,2000,50),mc.cores=14,
# 	FUN=function(x) { 
# 		kmodes(save.test, x, iter.max = 10000, weighted = FALSE, fast=T)
# 	})

stop()


# try PAM / medioids (what MetaBAT uses)

# traverse the tree to build bins by  merging?
# use completion as metric
# for manual clustering, pick a point, it picks the most common neighbors
# that cluster with that point, colored by frequency of co-occurrence in kmodes

# find my fastdist script

# can just do 100x kmeans3d and take the mode?

# daisy.out <- daisy(as.data.frame(save.class)[1:10000,],metric="gower")	 #slow, giant matrix
# pd.out <- parDist(as.matrix(daisy.out),method="binary")
# 
# pd.pca <- FastPCA(as.matrix(pd.out),50)
# pd.pca.uniq <- pd.pca$x
# plot(pd.pca.uniq[,1:2])


# rock fails
#rock.out <-rockCluster(as.matrix(save.test)[1:100,], 3)
# 
# for(j in 1:100) {
# 	for(i in 1:length(kmodes.out)) {
# 	try({
# 		plot3d(cov.tsne3d[which(kmodes.out[[i]]$cluster==j),],pch='.',bg=NULL,col=alpha("grey",0.4),size=0.2,main=count,xlim=range(cov.tsne3d[,1]),ylim=range(cov.tsne3d[,2]),zlim=range(cov.tsne3d[,3]),axes=F)
# 		plot3d(cov.tsne3d[which(kmodes.out[[i]]$cluster==j),],pch=21,col=alpha(kmodes.out[[i]]$cluster[which(kmodes.out[[i]]$cluster==j)],0.5),cex=4,main=count,add=T)
# 	})
# 		Sys.sleep(2)
# 	}
# }
# 
# #tsne 2d
# cov.tsne.list <- Rtsne(cov.pca.uniq, dims=2, max_iter=2000, verbose=TRUE, theta=0.1, pca=FALSE)
# cov.tsne <- cov.tsne.list$Y
# rownames(cov.tsne) <- rownames(covtab)



# 		print(paste0(count," ",m))
# 		tmp.class <- matrix(0,nrow=nrow(covtab),ncol=11)
# 		rownames(tmp.class) <- rownames(covtab)
# 		picked <- (fkm.tsne3d$class==m)
# 		plot3d(cov.tsne3d[picked,],pch=21,bg=NULL,col=alpha(fkm.tsne3d$class[picked],1),size=5,main=dms,add=T,xlim=range(cov.tsne3d[,1]),ylim=range(cov.tsne3d[,2]),zlim=range(cov.tsne3d[,3]))

# 		tmp.pca <- FastPCA(covtab[picked,],50)
# 		tmp.pca.uniq <- tmp.pca$x

# 		no need
#		tmp.tsne.list <- Rtsne(tmp.pca.uniq, dims=2, max_iter=1000, verbose=TRUE, theta=0.5, pca=FALSE)
#		tmp.tsne <- tmp.tsne.list$Y
#		rownames(tmp.tsne) <- rownames(covtab)[picked]

#		picked.len <- foam.meta[rownames(tmp.pca.uniq),"contigLen"]

# 		tmp.tsne3d.list <- Rtsne(tmp.pca.uniq, dims=3, max_iter=1000, verbose=F, theta=0.3, pca=FALSE)
# 		tmp.tsne3d <- tmp.tsne3d.list$Y
# 		rownames(tmp.tsne3d) <- rownames(covtab)[picked]

#clusgap fail
#		clusgap <- clusGap(cov.tsne3d[picked,], kmeans, 30, B = 100, verbose = interactive())

#	mcl fail
# 		tmp.mcl <- mcl(as.matrix(dist(cov.tsne3d[picked,])),addLoops=F,allow1=T)

# alas Mclust freezes up sometimes
# 		d_clust <- Mclust(as.matrix(cov.tsne3d[picked,]), G=1:40,verbose=F)
# # 		plot(d_clust,"BIC")
# 		print(paste0(count," ",m," ",dim(d_clust$z)))
# 		d_clust$classification <- d_clust$classification*m*dms
# 		return(d_clust$classification)

#	 	plot3d(cov.tsne3d,pch='.',add=T,bg=NULL,col=alpha("grey",0.4),size=0.2,main=count,xlim=range(cov.tsne3d[,1]),ylim=range(cov.tsne3d[,2]),zlim=range(cov.tsne3d[,3]),axes=F)
# 		plot3d(cov.tsne3d[picked,],pch=21,bg=NULL,col=alpha(d_clust$classification,0.5),size=5,main=dms,add=T,xlim=range(cov.tsne3d[,1]),ylim=range(cov.tsne3d[,2]),zlim=range(cov.tsne3d[,3]))

#		plot(tmp.tsne3d[,1:2],pch=21,col=NULL,bg=alpha("grey",1),cex=rescale(sqrt(picked.len),c(0.5,4)),main=paste0(m," md=",round(median(picked.len))," mn=",round(mean(picked.len))," N=",sum(picked)))
#		points(tmp.tsne3d[,1],tmp.tsne3d[,2],pch=21,col=NULL,bg=alpha(fkm.class,0.5),cex=rescale(sqrt(picked.len),c(0.5,4)),main=m)
#		points(tmp.tsne3d[,1],tmp.tsne3d[,2],pch=21,col=NULL,bg=alpha(d_clust$classification,0.5),cex=rescale(sqrt(picked.len),c(0.5,4)),main=m)

# 		loop <- 0
# 		for(l in 5:15) {
# 			loop <- loop+1
# 			count <- count+1
# 			fkm.tmp <- fastKmeans(tmp.tsne3d, l, iter.max = 10000, project = TRUE, threads = 0)
# 			fkm.class <- fkm.tmp$class
# 			rownames(fkm.class) <- rownames(tmp.tsne3d)
# 			tmp.class[rownames(fkm.class),loop] <- fkm.class
# 			rn <- rownames(fkm.class)
# 			print(l)
# 			for(cls in unique(fkm.class)) {
# 				whc <- which(fkm.class==cls)
# 				ln <- length(whc)
# # 				print(cls)
# 				for(i in 1:ln) {
# 					indi <- rn[whc[i]]
# 					for(j in i:ln) {
# 						indj <- rn[whc[j]]
# 						save.mat[indi,indj] <- save.mat[indi,indj]+1
# 					}
# 				}
# 			}
# 		}
# 		par(mfrow=c(1,2))
# 		plot(Rtsne(unique(FastPCA(save.mat,50)$x),pca=F)$Y)



# 	tmp.a <- as.data.frame(as.factor((do.call(c, save.tmp))))
# 	tmp.daisy <- daisy(as.data.frame(tmp.a),metric="gower")	 #slow, giant matrix

# 	graphics.off()
	


stop()

# 		nb <- NbClust(data=tmp.pca.uniq, diss=NULL, distance = "euclidean", 
# 				min.nc=2, max.nc=15, method = "kmeans", 
# 				index = "alllong", alphaBeale = 0.1)
# 		hist(nb$Best.nc[1,], breaks = max(na.omit(nb$Best.nc[1,])))


# bh-tsne 2D
cov.tsne.list <- Rtsne(cov.pca.uniq, dims=2, max_iter=2000, verbose=TRUE, theta=0.5, pca=FALSE)
saveRDS(cov.tsne.list,file="annotation/bins/tsne-embed2d.rds")
# cov.tsne.list <- readRDS(file="tsne-embed2d.rds")
cov.tsne <- cov.tsne.list$Y
rownames(cov.tsne) <- rownames(covtab)


fkm.tsne <- fastKmeans(cov.tsne, dims, iter.max = 1000, project = TRUE, threads = 0)
#plot(cov.tsne,pch=21,bg=NULL,col=alpha(fkm.tsne$class,0.4),cex=0.5,main=dims)
plot(cov.tsne,pch=21,col=NULL,bg=alpha(topo.colors(256)[round(rescale(cov.tsne3d[,1],c(1,256)))],0.4),cex=0.5)
plot(cov.tsne,pch=21,col=NULL,bg=alpha(fkm.tsne3d$class,0.4),cex=0.5,main=dims)

#UMAP dimensional reduction
#2d
cov.umap.list <- umap(cov.pca.uniq)
cov.umap <- cov.umap.list$layout
rownames(cov.umap) <- rownames(covtab)

plot(cov.umap,pch=21,col=NULL,bg=alpha(topo.colors(256)[round(rescale(cov.tsne3d[,1],c(1,256)))],0.4),cex=0.5,main=dims)

fkm.umap <- fastKmeans(cov.umap, dims, iter.max = 1000, project = TRUE, threads = 0)
plot3d(cov.umap,pch=21,bg=NULL,col=alpha(fkm.umap$class,0.4),cex=0.5,main=dims)

#3d
umap.configs <- umap.defaults
umap.configs$n_components=3
cov.umap3d.list <- umap(cov.pca.uniq, config=umap.configs)
cov.umap3d <- cov.umap3d.list$layout
rownames(cov.umap3d) <- rownames(covtab)

fkm.umap3d <- fastKmeans(cov.umap3d, 3, iter.max = 1000, project = TRUE, threads = 0)
plot3d(cov.umap3d,pch=21,bg=NULL,col=alpha(fkm.umap3d$class,0.4),cex=0.5,main=dims)


# binning
dims=100
km <- kmeans(cov.tsne,dims,nstart=20,iter.max=1000);

#plot(cov.tsne,pch=21,col=NULL,bg=alpha(topo.colors(256)[round(rescale(cov.umap[,3],c(1,256)))],0.4),cex=0.5,main="umap3")
#plot(cov.tsne,pch=21,col=NULL,bg=alpha(topo.colors(256)[round(rescale(cov.umap[,2],c(1,256)))],0.4),cex=0.5,main="umap2")
#plot(cov.tsne,pch=21,col=NULL,bg=alpha(topo.colors(256)[round(rescale(cov.umap[,1],c(1,256)))],0.4),cex=0.5,main="umap1")
par(mfrow=c(2,2))
for(i in 1:ncol(save.test)) {
	plot(cov.tsne,pch=21,col=NULL,bg=alpha(topo.colors(256)[round(rescale(cov.tsne3d[,3],c(1,256)))],0.4),cex=rescale(sqrt(foam.meta$contigLen),c(0.3,4)),main="tsne3")
	plot(cov.tsne,pch=21,col=NULL,bg=alpha(topo.colors(256)[round(rescale(cov.tsne3d[,2],c(1,256)))],0.4),cex=rescale(sqrt(foam.meta$contigLen),c(0.3,4)),main="tsne2")
	plot(cov.tsne,pch=21,col=NULL,bg=alpha(topo.colors(256)[round(rescale(cov.tsne3d[,1],c(1,256)))],0.4),cex=rescale(sqrt(foam.meta$contigLen),c(0.3,4)),main="tsne1")
	plot(cov.tsne,pch=21,col=NULL,bg=alpha(rainbow(length(unique(save.test[,i])))[as.factor(save.test[,i])],0.4),cex=rescale(sqrt(foam.meta$contigLen),c(0.3,4)),main="tsne1")
}


pdf(file="annotation/bins/cluster.pdf",width=24,height=24)
for(i in 1:ncol(save.test)) {
	plot(cov.tsne,pch=21,col=NULL,bg=alpha(rainbow(length(unique(save.test[,i])))[as.factor(save.test[,i])],0.4),cex=rescale(sqrt(foam.meta[rownames(covtab),"contigLen"]),c(0.1,5)),main="tsne1")
}
dev.off()
rownames(cov.tsne) <- rownames(covtab)


pdf(file="annotation/bins/cluster2.pdf",width=24,height=24)
siz <- rescale(sqrt(dts[rownames(cov.tsne),"V6"]),c(0.1,5))

for(i in 2:8) {
colr <- alpha(as.numeric(as.factor(gcut(dts[rownames(cov.tsne),"V2"]," ",i))),0.3)
plot(cov.tsne,pch=21,col=NULL,bg=colr, cex=siz,main=i)
}

colr <- alpha(as.numeric(as.factor(gcut(dts[rownames(cov.tsne),"V4"],".",4))),0.3)
plot(cov.tsne,pch=21,col=NULL,bg=colr, cex=siz,main=i)

dev.off()


# manual binning
tsne8 <- readRDS("tsne8.rds")
df.meta <- readRDS("df.meta.rds")
xmin <- min(tsne8$Y[,1])
xmax <- max(tsne8$Y[,1])
xwid <- (xmax - xmin)/1
ymin <- min(tsne8$Y[,2])
ymax <- max(tsne8$Y[,2])
ywid <- (ymax-ymin)/1

clength <- df.meta$contigLen[match(rownames(tsne8$Y),df.meta$node)]
plot(tsne8$Y,pch='.')

selectedPoints <- list()
i=1
for(x in 0:1) {
	for(y in 0:1) {
		for(i in 1:100) {
			plot(tsne8$Y, xlim=c(xmin + x*xwid,xmin + (x+1)*xwid), ylim=c(ymin+y*ywid, ymin + (y+1)*ywid), pch=21, col=NULL, bg=alpha("black",0.4), cex=rescale(sqrt(sqrt(clength)),to=c(0.1,5)))
			for(p in 1:length(selectedPoints)) {
				if(length(selectedPoints)<1) next
				bound <- attributes(selectedPoints[[p]])$gate
				lines(rbind(bound,bound[1,]))
			}
			try({selectedPoints[[(length(selectedPoints)+1)]] <- fhs(tsne8$Y, mark = TRUE)})
		}
	}
}


plot(cov.tsne,pch=21,col=NULL,bg=alpha(topo.colors(256)[round(rescale(cov.tsne3d[,1],c(1,256)))],0.4),cex=sqrt(foam.meta$contigLen)/1000,main="tsne1")

par(mfrow=c(3,7))
for(i in sort(unique(foam.meta$foam_1))) {
	plot(cov.tsne,pch=21,col=NULL,bg=alpha("grey",0.4),cex=sqrt(foam.meta$contigLen)/1000,main=i)
	picked1 <- foam.meta[rownames(cov.tsne),"foam_1"]==i
	picked2 <- foam.meta[rownames(cov.tsne),"contigLen"]>1
	picked3 <- picked1 & picked2
	points(cov.tsne[picked3,],pch=21,col=NULL,bg=alpha(fkm.tsne3d$class[picked3],0.4),cex=sqrt(foam.meta$contigLen[picked3])/500,main=i)
}

#foam.meta <- df.meta
foam.meta <- df.meta[!duplicated(df.meta$node),]
rownames(foam.meta) <- foam.meta$node

plot(tsne8$Y, xlim=c(xmin + x*xwid,xmin + (x+1)*xwid), ylim=c(ymin+y*ywid, ymin + (y+1)*ywid), pch=21, col=NULL, bg=alpha("black",0.4), cex=rescale(sqrt(sqrt(clength)),to=c(0.1,5)))

for(p in 1:length(selectedPoints)) {
	print(p)
	print(sum(foam.meta[selectedPoints[[p]],"contigLen"]))
	print(mean(foam.meta[selectedPoints[[p]],"totalAvgDepth"]))
#	write.table(file=paste0("bin",p,".txt"),cbind(paste0("bin",p),foam.meta[selectedPoints[[p]],"contig"]),sep="\t",quote=F,col.names=F,row.names=F)
	bound <- attributes(selectedPoints[[p]])$gate
	lines(rbind(bound,bound[1,]))
	text(bound[1,],cex=0.5, labels=paste0("bin",p,":",round(sum(foam.meta[selectedPoints[[p]],"contigLen"])/1e6,2),"Mbp"), col="green")
#	foam.meta[selectedPoints[[p]],]
	Sys.sleep(0.5)
}


saveRDS(file="selectedPoints.rds",selectedPoints)






# igraph binning

# this is to parallelize the intersect() distance
library(RcppArmadillo)
# Use RcppXPtrUtils for simple usage of C++ external pointers
library(RcppXPtrUtils)

# compile user-defined function and return pointer (RcppArmadillo is used as dependency)
euclideanFuncPtr <- cppXPtr(
	"double customDist(const arma::mat &A, const arma::mat &B) {
	return sqrt(arma::accu(arma::square(A - B)));
	return arma::intersect(A,B).n_elem;
	}", depends = c("RcppArmadillo"))

intersectFuncPtr <- cppXPtr(
	"double customDist(const arma::mat &A, const arma::mat &B) {
	arma::mat C;
	C = arma::intersect(A,B);
	return C.n_elem;
	}", depends = c("RcppArmadillo"))


save.df <- as.data.frame(save.test,stringsAsFactors=T)
save.num <- as.data.frame(sapply(save.df,as.numeric))
rownames(save.num) <- rownames(save.df)

# gets the number of duplicated rows
save.agg <- aggregate(list(numdup=rep(1,nrow(save.num))), save.num, length)
save.uniq <- unique(save.num)
sA <- apply(save.agg[,1:ncol(save.uniq)],1,paste,collapse=' ')
sU <- apply(save.uniq,1,paste,collapse=' ')
save.agg <- save.agg[match(sU,sA),]
rownames(save.agg) <- rownames(save.uniq)

#oops have 2 separate datasets loaded
foam.mm <- foam.meta[rownames(save.num),]

sN <- apply(save.num,1,paste,collapse=' ')
save.size <- aggregate(foam.mm$contigLen,list(y=sN),sum)
save.size <- save.size[match(sU,save.size$y),]
rownames(save.size) <- rownames(save.uniq)
save.size <- as.data.frame(save.size[,2])

dim(save.uniq)

#19000 row matrix takes about 3 minutes on 14 threads
#4000 takes 3 seconds
save.dist <- as.matrix(parDist(as.matrix(save.uniq), method="custom", func = intersectFuncPtr))
rownames(save.dist) <- rownames(save.df)[!duplicated(save.df)]
colnames(save.dist) <- rownames(save.df)[!duplicated(save.df)]

#10s of millions of edges have very few hits, let's get rid of them
h <- hist(save.dist[save.dist!=0],breaks=1:ncol(save.uniq))
plot(log10(h$counts))

save.dist[save.dist < 2] <- 0

# library(ggnet)

g1 <- graph_from_adjacency_matrix(save.dist, mode="upper", weighted=T,add.colnames=NULL)
g1 <- set_vertex_attr(graph=g1, name="size", value=save.size[,1])

paths.in <- read.table("spades_all_2018-10-03-10-44-21.paths.tsv",sep="\t")
g2 <- graph_from_data_frame(paths.in, directed = FALSE, vertices = NULL)
V(g2)$name)

cov.tsne.list <- readRDS("tsne-embed2d.rds")
cov.tsne <- cov.tsne.list$Y
rownames(cov.tsne) <- rownames(covtab)
l <- cov.tsne[rownames(save.dist),]

#somethings still wrong

pdf(file="graph.pdf",width=48,height=48)
for(i in rev(seq(min(E(g1)$weight),max(E(g1)$weight)))) {
	print(i)
	g1ns <- igraph::delete.edges(g1, which(E(g1)$weight < i))
	plot(g1ns, layout=l, vertex.label=NA, edge.arrow.size=0, edge.weights=E(g1ns)$weight/50, vertex.size=sqrt(as.numeric(V(g1ns)$size))/300)
	plot(g2,layout=ll
}
dev.off()




# for(i in seq_along(Graph)) {
#   assign(paste0('g', i), Graph[i])
# }
# for(j in 1:ncol(x)) x[[j]][is.na(x[[j]])] = 0



#heatmap of aggregate bin relative abundance in samples


