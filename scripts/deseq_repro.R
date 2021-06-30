wd <- "~/Documents/TAAAAF/Stage/liver/"
image_dir <- file.path(wd, "analysis/R_images/")
setwd(wd)

## load required packages
library(BiocManager)
library(DESeq2)
library(RColorBrewer)
library(calibrate)
library(gplots)
library(genefilter)
#library(vsn) for microarray data, not sure if useful/necessary
library(pheatmap)
library(tximeta)
library(magrittr)
library(stringr)
library("AnnotationDbi")
library("org.Hs.eg.db")
library(ggvenn)


#################################
#   Importation                 #
#################################
dir_quants <- file.path(wd, "data/quants")
runTable <- read.table("SraRunTable.txt", header = T, sep = ",")

#Treatment is not really suitable as a condition, but how to automate "shrinkage"?
runTable$Treatment <- sapply(runTable$Treatment, str_trunc, width = 7, side = "right", ellipsis = "")
coldata <- data.frame(names = runTable$Run, condition = runTable$Treatment)
coldata$files <- file.path(dir_quants, coldata$names, "quant.sf")
file.exists(coldata$files)

#Annotation et importation
se <- tximeta(coldata)
gse <- summarizeToGene(se)

#Arrengement des facteurs
gse$condition %<>% factor()
gse$condition %<>% relevel("vehicle") #IMPORTANT
levels(gse$condition)[levels(gse$condition)=="recombi"] <- "PDGF_BB"
levels(gse$condition)[levels(gse$condition)=="TGF-bet"] <- "TGF_beta"

################################
#     Exploratory Analysis     #
################################

dds <- DESeqDataSet(gse, design = ~ condition) 


#Basic filter
#Only keeps genes with more than count across fragment
keep <- rowSums(counts(dds)) > 1 
dds <- dds[keep,]
nrow(dds)

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)

#Heatmap of Sample-Sample distance
condition <- dds$condition
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])
sampleDists <- as.matrix(dist(t(assay(rld))))
png(file.path(image_dir,"qc-heatmap-samples.png"), w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "red", "blue"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()



# Principal components analysis
#Soham's code : 
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL,
                     legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]),
                     1, paste, collapse = " : "))
  
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21,
       xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2,
                                    labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)             
}

png(file.path(image_dir,"qc-pca.png"), 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
dev.off()

#DESeq2 built-in :
#I prefer this one cause it uses the entire space without using a xlim entry
#So it is more portable
svg(file.path(image_dir,"DESeq2-pca.svg"), 7, 7, pointsize=20)
plotPCA(rld)
dev.off()





########################################
#           DEG ANALYSIS               #
########################################

dds <- DESeq(dds)

# Plot dispersions
png(file.path(image_dir, "qc-dispersions.png"), 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

#One should think about a way to automate the selection of contrast
res_TGF <- results(dds, contrast = c("condition", "TGF_beta", "vehicle"))
res_PDGF <- results(dds, contrast = c("condition", "PDGF_BB", "vehicle"))
res_P_T <- results(dds, contrast = c("condition", "PDGF_BB", "TGF_beta"))


#/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
#-------------------------------------------------------------------
#             I used GRCH38 and not hg19 as in the paper

#Annotation (not sure it works well)
ens.str <- substr(rownames(res_TGF), 1, 15)
res_TGF$symbol <- res_PDGF$symbol <- 
  res_P_T$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

#Finding which genes are upregulated vs vehicle in both treatment
table(res_TGF$padj < 0.05 & res_TGF$log2FoldChange > 1.5)
table(res_PDGF$padj < 0.05 & res_PDGF$log2FoldChange > 1.5)
sig_TGF <- which(res_TGF$padj < 0.05 & res_TGF$log2FoldChange > 1.5)
sig_PDGF <- which(res_PDGF$padj < 0.05 & res_PDGF$log2FoldChange > 1.5)
both_sig <- intersect(sig_PDGF, sig_TGF)

#drawing a Venn diagramm of upregulated genes
list_sig <- list("pdgf" = sig_PDGF, "tgf" = sig_TGF)
ggvenn(list_sig, c("pdgf", "tgf"))

#finding the names of co-up-regulated genes
both_sig_names <- res_TGF[both_sig, "symbol"]


#---------Compare with paper selection V1 ENSEMBL -> Symbol
genes_repro <- read.table(file.path(wd, "data/Repro_data/genes_up_propre.txt"))
genes_repro_vec <- genes_repro$V1
#liste des noms en commun entre les deux 
list_names_comp <- list("mine"=both_sig_names, "paper"=genes_repro_vec)
ggvenn(list_names_comp, c("mine", "paper"))
#not looking great are we, only 67% of 112 genes found.

#--------Compare with paper selection V2 Symbol -> ENSEMBL
#The following line assigns duplicate ENS ref 
genes_repro_ens <- mapIds(org.Hs.eg.db,
                         keys=unname(genes_repro_vec),
                         column="ENSEMBL",
                         keytype="SYMBOL",
                         multiVals="first")
both_sig_ens <- rownames(res_TGF)[both_sig]
list_names_comp_ens <- list("mine"=both_sig_ens, "paper"=unname(genes_repro_ens))

ggvenn(list_names_comp_ens, c("mine", "paper"))

###################################
#      Looking at indiv genes     #
###################################
#EZH2
lfc_tgf <- res_TGF[which(res_TGF$symbol=="EZH2"),]$log2FoldChange
lfc_pdgf <- res_PDGF[which(res_PDGF$symbol=="EZH2"),]$log2FoldChange
fc_tgf <- 2**lfc_tgf
fc_pdgf <- 2**lfc_pdgf
#Here I'm good, it's very close to the paper's numbers

pdgfa <- res_TGF[which(res_TGF$symbol=="PDGFA"),]$log2FoldChange
pdgfb <- res_TGF[which(res_TGF$symbol=="PDGFB"),]$log2FoldChange
vegfa <- res_TGF[which(res_TGF$symbol=="VEGFA"),]$log2FoldChange
fgf2 <- res_TGF[which(res_TGF$symbol=="FGF2"),]$log2FoldChange
#Close to the paper, too
ctgf <- res_TGF[which(res_TGF$symbol=="CCN2"),]$log2FoldChange
#needed to change name



### GO analysis
library("biomaRt")
library("goseq")
annotation <- read.delim("Annot.txt", header=TRUE)
diffexp <- read.delim("diffexpr-results.csv", sep = ",")

matched <- merge(annotation, diffexp, by.x="Gene_ID", by.y="Gene")
pwf <- nullp(matched, "Name", "Gene", bias.data = shrinkLvV$medianTxLength)





