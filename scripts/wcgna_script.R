#setwd("~/Documents/Collaborations/Mariana_Franck/______NEW ANALYSIS/WGCNA")

library(WGCNA)
library(flashClust)
library(gplots)
library(ggplot2)

allowWGCNAThreads()
#enableWGCNAThreads()
options(stringsAsFactors = FALSE)
countdata <- read.delim("diffexpr-results_annot.csv", sep = ",")
head(countdata)

data1 <- countdata[,c(12:25, 27:29)]
data1$variance = apply(data1, 1, var)
index = data1$variance >= quantile(data1$variance, c(.90))
index[17501] = TRUE
data2 = data1[data1$variance >= quantile(data1$variance, c(.90)), ] #10% most variable genes
data2$variance <- NULL 
data3 <- t(data2)
data4<- apply(data3, 2, as.numeric)


dat = as.data.frame(t(countdata[, c(2)]))
datExpr = dat[index]
#View(datExpr)
head(datExpr)
dim(datExpr)

powers = c(c(1:20), seq(from = 12, to=20, by=2))# in practice this should include powers up to 20.
sft<- pickSoftThreshold(data4,dataIsExpr = TRUE,
                        powerVector = powers,corFnc = cor,
                        corOptions = list(use = 'p'),networkType = "signed hybrid")


#plot SFT
sizeGrWindow(9, 7)
par(mfrow = c(1,2))
cex1 = 0.9
png("___softpower.png",width =10,height=10,unit="cm",res=100)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2", main = paste("Scale independence"))
abline(h=0.9,col="red")# Red line corresponds to using an R^2 cut-off
dev.off()
png("___scale connectivity.png",width =10,height=10,unit="cm",res=100)
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

softPower = 13
adj= adjacency(data4,type = "signed hybrid", power = softPower)
#turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(data4,networkType = "signed hybrid", TOMType = "signed", power = softPower)
colnames(TOM) =rownames(TOM) = datExpr
dissTOM=1-TOM
#Hierarchical clustering of the genes based on the TOM dissimilarity measure
geneTree = flashClust(as.dist(dissTOM),method="average")
#Plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="",cex=0.0013)


minModuleSize <- 25
ds <- 3
cutHeight <- 0.99999
dthresh <- 0.12 

# Module identification using dynamic tree cut
dynamicMods = cutreeHybrid(dendro = geneTree,cutHeight = cutHeight,distM = dissTOM,
                           deepSplit = ds, pamRespectsDendro = TRUE,pamStage = TRUE, 
                           minClusterSize = minModuleSize);
table(dynamicMods[[1]])
merged <- mergeCloseModules(exprData = data4,colors =dynamicMods$labels ,
                            cutHeight = dthresh)
summary(merged)

mColorh <- labels2colors(merged$colors)
mLabelh <- "DS=2,MMS=200,DCOR=0.15"

table(mColorh)
plotDendroAndColors(dendro = geneTree, mColorh, groupLabels = mLabelh,
                    addGuide=TRUE,dendroLabels=FALSE,main="Dendrogram With Defined Modules")
mergedColors = merged$colors
mergedMEs = t(merged$newMEs)


moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
#table(tree1)
dynamicColors = labels2colors(merged$colors)
table(dynamicColors)

png("___dendrogram.png",width =10,height=10,unit="cm",res=100)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, 
                    guideHang = 0.2, main = "Gene dendrogram and module colors")
dev.off()

restGenes= (dynamicColors != "grey")
diss1=1-TOMsimilarityFromExpr(data4[,restGenes], power = softPower)
collectGarbage()

colnames(diss1) =rownames(diss1) = datExpr[restGenes]
hier1=flashClust(as.dist(diss1), method="average" )
plotDendroAndColors(hier1, dynamicColors[restGenes], 
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, 
                    guideHang = 0.2, main = "HCLUST_Gene dendrogram and module colors")

diag(diss1) = NA
sizeGrWindow(7,7)
collectGarbage()

png("___TOMplot_90.png",width =10,height=10,unit="cm",res=100)
TOMplot(diss1, hier1, as.character(dynamicColors[restGenes]))
dev.off()
module_colors= setdiff(unique(dynamicColors), "grey")
for (color in module_colors){
  module=data4[which(dynamicColors==color)]
  module1 = datExpr[dynamicColors==color]
  write.table(t(module), paste("tmod_",color, ".txt",sep=""), sep="\t", 
              row.names=TRUE, col.names=TRUE,quote=FALSE)
}

module.order <- unlist(tapply(1:ncol(data4),as.factor(dynamicColors),I))
x<-list(module.order)

write.table(t(x), paste("x_module_90", ".txt",sep=""), sep="\t", 
            row.names=TRUE, col.names=TRUE,quote=FALSE)


m<-t(data4[,module.order])/apply(data4[,module.order],2,max)
names1 <- datExpr[,1]
colnames(m) <- names1
heatmap(m,zlim=c(0,1),col= redWhiteGreen(256),Rowv=NA,Colv=NA,labRow=NA,scale="row",RowSideColors=dynamicColors[module.order])
dev.off()

library(ComplexHeatmap)
library(circlize)

png("data_90.png",width =10,height=10,unit="cm",res=100)
p2  <-   Heatmap(m, name = "Row Z score"
                 ,col=colorRamp2(c(-1,-0.5,0,0.2,0.4,0.6,1), c("yellow","#fa5f5f","white","darkblue","blue","blue","#ffd14c"))
                 #,col=colorRamp2(c(-1,0,1),c("blue","white","red"))
                 ,cluster_rows = FALSE
                 ,clustering_distance_rows = "euclidean"
                 ,clustering_method_rows = "complete"
                 ,cluster_columns = FALSE
                 ,clustering_distance_columns = "euclidean"
                 ,clustering_method_columns = "complete"
                 #,rect_gp = gpar(col = "black"
                 ,na_col = "white"
                 ,km =1
                 ,heatmap_legend_param = list(at=c(0.6,1),color_bar = "continuous")
                 ,show_row_names = FALSE
                 ,show_row_dend = FALSE
                 #,column_names_gp = gpar(cex=0.1, font=0.1, col= "green"),
                 ,row_names_gp =  gpar(fontsize=12,fontcolor= "red"),
                 column_names_gp = gpar(fontsize = 20, fontface="bold")
                 ,show_column_names = FALSE)

module_color <- as.data.frame(dynamicColors[module.order])
colnames(module_color) <- "color"
c <- list(color=unique(module_color$color))
names(c$color) <- unique(module_color$color)
ha_row = rowAnnotation(df = module_color, col= c,width = unit(1, "cm"))
p2+ha_row

dev.off()

m1<-t(data4[,module.order])/apply(data4[,module.order],2,max)

ids = countdata[index, 1]
ID_corr =ids[as.numeric(rownames(m1))]
module_corrected<-data.frame(ID_corr, m)

ids2 = countdata[index, 2:11]
module_corrected1<-data.frame(ids2, module_corrected)

expt = data.frame(colnames(countdata[c(12:29)]))
expt = data.frame("genes", t(expt))
module_corrected<-data.frame(ID_corr, m)

write.table(module_corrected1, paste("___zscores_90_ALL", ".txt",sep=""), sep="\t", 
            row.names=TRUE, col.names=TRUE,quote=TRUE)

write.table(x, paste("variance", ".txt"))


MEList = moduleEigengenes(datExpr, colors = dynamicColors,nPC = 10)
MEs = MEList$eigengenes
png("Eigennetwork.png",width =10,height=10,unit="cm",res=100)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))
signif(cor(MEs,use = "p"),2)
dev.off()

nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
datTraits <- countdata$Name
head(datTraits)

MEs0 = moduleEigengenes(datExpr, dynamicColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);



##Make TREE Network
library(igraph)

#prepare data all datapoints
# Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
g <- graph.adjacency(
  as.matrix(as.dist(cor(dissTOM, method="pearson"))),
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)
# Simplfy the adjacency object
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE)
# Colour negative correlation edges as blue
E(g)[which(E(g)$weight<0)]$color <- "darkblue"
# Colour positive correlation edges as red
E(g)[which(E(g)$weight>0)]$color <- "darkred"
# Convert edge weights to absolute values
E(g)$weight <- abs(E(g)$weight)
# Remove edges below absolute Pearson correlation 0.8
g <- delete_edges(g, E(g)[which(E(g)$weight<0.8)])
# Remove any vertices remaining that have no edges
g <- delete.vertices(g, degree(g)==0)
# Assign names to the graph vertices (optional)
V(g)$name <- V(g)$name
# Change shape of graph vertices
V(g)$shape <- "sphere"
# Change colour of graph vertices
V(g)$color <- "skyblue"
# Change colour of vertex frames
V(g)$vertex.frame.color <- "white"
# Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
# Multiply scaled vales by a factor of 10
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
vSizes <- (scale01(apply(mergedMEs, 1, mean)) + 1.0)*2
# Amplify or decrease the width of the edges
edgeweights <- E(g)$weight * 2.0

# Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst <- mst(g, algorithm="prim")

# Plot the tree object
png("___network_all_99_diss1.png",width =50,height=50,unit="cm",res=300)
plot(
  mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.1,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="Network")
dev.off()

#mst.communities <- cluster_fast_greedy(mst)
mst.communities <- edge.betweenness.community(mst, weights=NULL, directed=FALSE)
mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
V(mst)$color <- mst.communities$membership + 1

png("___network_annot_all_99.png",width =10,height=10,unit="cm",res=600)
plot(
  mst.clustering, mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes/4,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.00001,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="Network")
dev.off()



#prepare data all cluster datapoints
# Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
g <- graph.adjacency(
  as.matrix(as.dist(cor(t(mergedMEs), method="pearson"))),
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)
# Simplfy the adjacency object
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE)
# Colour negative correlation edges as blue
E(g)[which(E(g)$weight<0)]$color <- "darkblue"
# Colour positive correlation edges as red
E(g)[which(E(g)$weight>0)]$color <- "darkred"
# Convert edge weights to absolute values
E(g)$weight <- abs(E(g)$weight)
# Remove edges below absolute Pearson correlation 0.8
g <- delete_edges(g, E(g)[which(E(g)$weight<0.5)])
# Remove any vertices remaining that have no edges
g <- delete.vertices(g, degree(g)==0)
# Assign names to the graph vertices (optional)
V(g)$name <- V(g)$name
# Change shape of graph vertices
V(g)$shape <- "sphere"
# Change colour of graph vertices
V(g)$color <- "skyblue"
# Change colour of vertex frames
V(g)$vertex.frame.color <- "white"
# Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
# Multiply scaled vales by a factor of 10
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
vSizes <- (scale01(apply(mergedMEs, 1, mean)) + 1.0) * 10

# Amplify or decrease the width of the edges
edgeweights <- E(g)$weight * 2.0
# Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst <- mst(g, algorithm="prim")
# Plot the tree object
png("___network_cluster.png",width =50,height=50,unit="cm",res=300)
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
  main="Network")
dev.off()

mst.communities <- cluster_fast_greedy(mst)
#mst.communities <- edge.betweenness.community(mst, weights=NULL, directed=FALSE)
mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
V(mst)$color <- mst.communities$membership + 1

png("___network_annot_cluster.png",width =10,height=10,unit="cm",res=300)
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
  main="Network")
dev.off()