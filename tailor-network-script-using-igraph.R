# These arguments are passed in on the command line via launchR.eachDiffDir.sh
# args:
# diffDir=\"${diffDir}\" inDir=\"${INPUT}/${diffDir}\" outDir=\"${OUTPUT}/${diffDir}\"
# under=\"${under}\" over=\"${over}\" Rplots=\"$\{OUTPUT}/${diffDir}/${R_PLOTS}\"
# FPKMmatrix=\"${OUTPUT}/${diffDir}/${FPKM_MATRIX}\" DiffTable=\"${OUTPUT}/${diffDir}/${DIFF_TABLE}\"
# under="CTRL";over="LUTS"

############################################################################
### code chunk number 1: Load Libraries and initialize cuff object
############################################################################
library(cummeRbund);library(knitr);library(igraph)
library(ALL);library(rpart);library(rpart.plot)
library(ape);library(lmtest);library(multtest)
library(randomForest);library(outliers)
library(msigdbr)
cuff=readCufflinks()
cuff
runInfo(cuff)
replicates(cuff)
replicates(cuff)$rep_name
over<-gsub("_[0:9]", "", replicates(cuff)$rep_name[1])
under<-gsub("_[0:9]", "", replicates(cuff)$rep_name[(length(replicates(cuff)$rep_name)/2)+1])
over;under
###################################################################################################
### code chunk number : Step 2, Read in, normalise, and identify genes with significant effects
####################################################################################################

# Sometimes diffData() doesnt work at gene level and the analysis must be performed at isoform level
# This requires mapping isoform IDs to gene names

gene_diff_data<-diffData(genes(cuff))
head(gene_diff_data)
g.rep.matrix<-repFpkmMatrix(genes(cuff))
head(g.rep.matrix)

#fpkm.file = file.path(FPKMmatrix)
rep.fpkm.file = file.path("repFPKMmatrix")
write.table(g.rep.matrix,file = rep.fpkm.file, sep = "  ", row.names = T, col.names = T,quote = F)

o.group<-grep(pattern = over,colnames(g.rep.matrix),ignore.case = T)
u.group<-grep(pattern = under,colnames(g.rep.matrix),ignore.case = T)

sig_gene_data<-subset(gene_diff_data, (significant =='yes'))
head(sig_gene_data)
nrow(sig_gene_data)

sig.g.rep.matrix<-g.rep.matrix[which(row.names(g.rep.matrix) %in% sig_gene_data$gene_id),]
head(sig.g.rep.matrix)
dim(sig.g.rep.matrix)

sig_gene_data<-subset(gene_diff_data, log2_fold_change > 10 | log2_fold_change < -10)
head(sig_gene_data)
nrow(sig_gene_data)


sig_over_gene_data<-subset(gene_diff_data, (significant =='yes') & (log2_fold_change > 0))
head(sig_over_gene_data)
nrow(sig_over_gene_data)

sig_under_gene_data<-subset(gene_diff_data, (significant =='yes') & (log2_fold_change < 0))
head(sig_under_gene_data)
nrow(sig_under_gene_data)

sig_over_gene_data<-subset(gene_diff_data, log2_fold_change > 10)
head(sig_over_gene_data)
nrow(sig_over_gene_data)

sig_under_gene_data<-subset(gene_diff_data, log2_fold_change < -10)
head(sig_under_gene_data)
nrow(sig_under_gene_data)

sig.g.o.rep.ma<-g.rep.matrix[which(row.names(g.rep.matrix) %in% sig_over_gene_data$gene_id),]
head(sig.g.o.rep.ma)
dim(sig.g.o.rep.ma)

sig.g.u.rep.ma<-g.rep.matrix[which(row.names(g.rep.matrix) %in% sig_under_gene_data$gene_id),]
head(sig.g.u.rep.ma)
dim(sig.g.u.rep.ma)


########################################################################
# If the code that was just run failed
# do this
# otherwise skip this
########################################################################

#iso_diff_data<-diffData(isoforms(cuff))
#i.rep.matrix<-repFpkmMatrix(isoforms(cuff))
#gene.xloc.matrix<-featureNames(isoforms(cuff))
#gene.list<-getGenes(cuff,geneId = gene.xloc.matrix)
#gene_annotation_data<-featureNames(gene.list@isoforms)
#gene_annotation_data<-gene_annotation_data[order(gene_annotation_data$tracking_id),]
g.rep.matrix<-repFpkmMatrix(genes(cuff))
gene.xloc.matrix<-featureNames(genes(cuff))
#gene.xloc.matrix<-featureNames(isoforms(cuff))
gene.list<-getGenes(cuff,geneId = gene.xloc.matrix)
gene_annotation_data<-featureNames(gene.list)
gene_annotation_data<-gene_annotation_data[order(gene_annotation_data$tracking_id),]
gene_rep_matrix<-cbind2(gene_annotation_data, g.rep.matrix)
gene_mapping=gene_annotation_data[,"gene_short_name"]
names(gene_mapping)=gene_annotation_data[,"tracking_id"]
ids=unique(gene_rep_matrix[,"tracking_id"])
gene_names=gene_mapping[ids]
gene.ids=unique(gene_rep_matrix[,"gene_short_name"])
gene_names=gene_mapping[gene.ids]
gene.rep.matrix=gene_rep_matrix[gene.ids,c(-1,-2)]
row.names(gene.rep.matrix)=gene.ids
sig.isos<-which(gene_annotation_data$tracking_id %in% sig_gene_data$isoform_id)
sig_gene_annotation_data<-gene_annotation_data[sig.isos,]
sig.rep.fpkm.file = file.path("sig.repFPKMmatrix")
write.table(sig.g.rep.matrix,file = sig.rep.fpkm.file, sep = "  ", row.names = T, col.names = T,quote = F)
over.group<-grep(pattern = over,colnames(sig.g.rep.matrix),ignore.case = T)
under.group<-grep(pattern = under,colnames(sig.g.rep.matrix),ignore.case = T)
rep.gene.counts<-repCountMatrix(genes(cuff))
sig.rep.gene.counts<-rep.gene.counts[which(row.names(rep.gene.counts) %in% sig_gene_data$gene_id),]
rep.count.file = file.path("repCOUNTmatrix")
write.table(rep.gene.counts,file = rep.count.file, sep = "  ", row.names = T, col.names = T,quote = F)
sig.rep.count.file = file.path("sig.repCOUNTmatrix")
write.table(sig.rep.gene.counts,file = sig.rep.count.file, sep = "  ", row.names = T, col.names = T,quote = F)

sig_over_isos<-which(gene_annotation_data$tracking_id %in% sig_over_gene_data$isoform_id)
sig_over_annotation_data<-gene_annotation_data[sig_over_isos,]

sig_under_isos<-which(gene_annotation_data$tracking_id %in% sig_under_gene_data$isoform_id)
sig_under_annotation_data<-gene_annotation_data[sig_under_isos,]

sig.g.o.rep.matrix<-g.rep.matrix[which(row.names(g.rep.matrix) %in% sig_over_annotation_data$gene_short_name),]
sig.o.rep.fpkm.file = file.path("sig.over.repFPKMmatrix")
write.table(sig.g.o.rep.matrix,file = sig.o.rep.fpkm.file, sep = "  ", row.names = T, col.names = T,quote = F)
over.g<-grep(pattern = over,colnames(sig.g.o.rep.matrix),ignore.case = T)
sig.g.o.rep.ma<-sig.g.o.rep.matrix[over.g]
sig.g.u.rep.matrix<-g.rep.matrix[which(row.names(g.rep.matrix) %in% sig_under_annotation_data$gene_short_name),]
sig.u.rep.fpkm.file = file.path("sig.under.repFPKMmatrix")
write.table(sig.g.u.rep.matrix,file = sig.u.rep.fpkm.file, sep = "  ", row.names = T, col.names = T,quote = F)
under.g<-grep(pattern = under,colnames(sig.g.u.rep.matrix),ignore.case = T)
sig.g.u.rep.ma<-sig.g.u.rep.matrix[under.g]

########################################################################
### code chunk number 4: Write adj mtrix and plot graph for over genes
########################################################################
network.file = file.path("Gene_Expression_Network_Analysis.pdf")
pdf(file=network.file)
#sig.g.rep.matrix
#Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
graph <- graph.adjacency(as.matrix(as.dist(cor(t(sig.g.rep.matrix), method="pearson"))),mode="undirected", weighted=TRUE, diag=FALSE)
head(graph)
#NB - Euclidean distances can also be used instead of correlation, but this changes the tutorial slightly:
#euclidean
#graph <- graph.adjacency(as.matrix(dist(sig.g.rep.matrix, method="euclidean")),
#                      mode="undirected", weighted=TRUE, diag=FALSE)

#Simplfy the adjacency object
graph <- simplify(graph, remove.multiple=TRUE, remove.loops=TRUE)
#Colour negative correlation edges as blue
E(graph)[which(E(graph)$weight<0)]$color <- "darkblue"
#Colour positive correlation edges as red
E(graph)[which(E(graph)$weight>0)]$color <- "darkred"
#Convert edge weights to absolute values
E(graph)$weight <- abs(E(graph)$weight)
#Change arrow size
#For directed graphs only
#E(over.graph)$arrow.size <- 1.0
#Remove edges below absolute Pearson correlation 0.8
graph <- delete_edges(graph, E(graph)[which(E(graph)$weight<0.8)])
#Assign names to the graph vertices (optional)
V(graph)$name <- V(graph)$name
#Change shape of graph vertices
V(graph)$shape <- "sphere"
#Change colour of graph vertices
V(graph)$color <- "skyblue"
#Change colour of vertex frames
V(graph)$vertex.frame.color <- "white"

#Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
#Multiple scaled vales by a factor of 10
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
vSizes <- (scale01(apply(sig.g.rep.matrix, 1, mean)) + 1.0) * 10

#Amplify or decrease the width of the edges
edgeweights <- E(graph)$weight * 2.0
#Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst <- mst(graph, algorithm="prim")
mst <- mst(graph)
#Plot the tree object
plot(mst, layout=layout.fruchterman.reingold, edge.curved=TRUE,
     vertex.size=vSizes, vertex.label.dist=-0.5, vertex.label.color="black",
     asp=FALSE, vertex.label.cex=0.6, edge.width=edgeweights, edge.arrow.mode=0,
     main="All correlated highly expressed genes included in network")

#################################################################################################
### code chunk number 5: Identify communities in the tree object based on 'edge betweenness'
###############################################################################################

mst.communities <- edge.betweenness.community(mst, directed=FALSE)
mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
V(mst)$color <- mst.communities$membership + 1

plot(
  mst.clustering, mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="All correlated highly expressed genes included in network")
)

plot(
  mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="All correlated highly expressed genes included in network")
)


########################################################################
### code chunk number 5: Write adj mtrix and plot graph for over genes
########################################################################

#Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
over.graph <- graph.adjacency(as.matrix(as.dist(cor(t(sig.g.o.rep.ma), method="pearson"))),
                              mode="undirected", weighted=TRUE, diag=FALSE)

head(over.graph)
#NB - Euclidean distances can also be used instead of correlation, but this changes the tutorial slightly:

#over.graph <- graph.adjacency(as.matrix(dist(sig.g.o.rep.ma, method="euclidean")),
#                      mode="undirected", weighted=TRUE, diag=FALSE)
#Simplfy the adjacency object
over.graph <- simplify(over.graph, remove.multiple=TRUE, remove.loops=TRUE)
#Colour negative correlation edges as blue
E(over.graph)[which(E(over.graph)$weight<0)]$color <- "darkblue"
#Colour positive correlation edges as red
E(over.graph)[which(E(over.graph)$weight>0)]$color <- "darkred"
#Convert edge weights to absolute values
E(over.graph)$weight <- abs(E(over.graph)$weight)
#Change arrow size
#For directed graphs only
#E(over.graph)$arrow.size <- 1.0

#Remove edges below absolute Pearson correlation 0.8
over.graph <- delete_edges(over.graph, E(over.graph)[which(E(over.graph)$weight<0.8)])
#Assign names to the graph vertices (optional)
V(over.graph)$name <- V(over.graph)$name
#Change shape of graph vertices
V(over.graph)$shape <- "sphere"
#Change colour of graph vertices
V(over.graph)$color <- "skyblue"
#Change colour of vertex frames
V(over.graph)$vertex.frame.color <- "white"

#Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
#Multiple scaled vales by a factor of 10
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
vSizes <- (scale01(apply(sig.g.o.rep.ma, 1, mean)) + 1.0) * 10

#Amplify or decrease the width of the edges
edgeweights <- E(over.graph)$weight * 2.0
#Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst <- mst(over.graph, algorithm="prim")

#Plot the tree object
plot(mst, layout=layout.fruchterman.reingold, edge.curved=TRUE,
     vertex.size=vSizes, vertex.label.dist=-0.5, vertex.label.color="black",
     asp=FALSE, vertex.label.cex=0.6, edge.width=edgeweights, edge.arrow.mode=0,
     main=substitute(paste(o," patient significant gene network"), list(o=eval(over))))

#################################################################################################
### code chunk number 6: Identify communities in the tree object based on 'edge betweenness'
###############################################################################################

mst.communities <- edge.betweenness.community(mst, directed=FALSE)
mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
V(mst)$color <- mst.communities$membership + 1

plot(
  mst.clustering, mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main=substitute(paste(o," patient significant gene network"), list(o=eval(over))))

plot(
  mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main=substitute(paste(o," patient significant gene network"), list(o=eval(over))))

########################################################################
### code chunk number 5: Write adj mtrix and plot graph for under genes
########################################################################

#Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
under.graph <- graph.adjacency(as.matrix(as.dist(cor(t(sig.g.u.rep.ma), method="pearson"))),
                              mode="undirected", weighted=TRUE, diag=FALSE)

#NB - Euclidean distances can also be used instead of correlation, but this changes the tutorial slightly:
#g <- graph.adjacency(as.matrix(dist(estrogenMainEffects, method="euclidean")),
#mode="undirected", weighted=TRUE, diag=FALSE)
#Simplfy the adjacency object
under.graph <- simplify(under.graph, remove.multiple=TRUE, remove.loops=TRUE)
#Colour negative correlation edges as blue
E(under.graph)[which(E(under.graph)$weight<0)]$color <- "darkblue"
#Colour positive correlation edges as red
E(under.graph)[which(E(under.graph)$weight>0)]$color <- "darkred"
#Convert edge weights to absolute values
E(under.graph)$weight <- abs(E(under.graph)$weight)
#Change arrow size
#For directed graphs only
#E(over.graph)$arrow.size <- 1.0

#Remove edges below absolute Pearson correlation 0.8
under.graph <- delete_edges(under.graph, E(under.graph)[which(E(under.graph)$weight<0.8)])
#Assign names to the graph vertices (optional)
V(under.graph)$name <- V(under.graph)$name
#Change shape of graph vertices
V(under.graph)$shape <- "sphere"
#Change colour of graph vertices
V(under.graph)$color <- "skyblue"
#Change colour of vertex frames
V(under.graph)$vertex.frame.color <- "white"

#Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
#Multiple scaled vales by a factor of 10
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
vSizes <- (scale01(apply(sig.g.u.rep.ma, 1, mean)) + 1.0) * 10

#Amplify or decrease the width of the edges
edgeweights <- E(under.graph)$weight * 2.0
#Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst <- mst(under.graph, algorithm="prim")

#Plot the tree object
plot(mst, layout=layout.fruchterman.reingold, edge.curved=TRUE,
     vertex.size=vSizes, vertex.label.dist=-0.5, vertex.label.color="black",
     asp=FALSE, vertex.label.cex=0.6, edge.width=edgeweights, edge.arrow.mode=0,
     main=substitute(paste(u," patient significant gene network"), list(u=eval(under))))


#################################################################################################
### code chunk number 4: Identify communities in the tree object based on 'edge betweenness'
###############################################################################################

mst.communities <- edge.betweenness.community(mst, directed=FALSE)
mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
V(mst)$color <- mst.communities$membership + 1


plot(mst.clustering, mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main=substitute(paste(u," patient significant gene network"), list(u=eval(under))))

plot(
  mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main=substitute(paste(u," patient significant gene network"), list(u=eval(under))))

#################################################################################################
### code chunk number 5: Step 4, Further analyses
###############################################################################################

#Check the vertex degree, i.e., number of connections to each vertex
degree(mst)

#Output information for each community, including vertex-to-community assignments and modularity
commSummary <- data.frame(
  mst.communities$names,
  mst.communities$membership,
  mst.communities$modularity
)
colnames(commSummary) <- c("Gene", "Community", "Modularity")
options(scipen=999)
commSummary

#Compare community structures using the variance of information (vi) metric (not relevant here)
#Community structures that are identical will have a vi=0
compare(mst.communities, mst.communities, method="vi")
dev.off()
