# load libraries -----------------------------------------------------------------
library(shiny)
library(BiocManager)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tximeta)
library(tximport)
library(magrittr)
library(stringr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(shinythemes)
library(GWENA)
library(gprofiler2)
library(ggplot2)
library(plotly)
library(shinyvalidate)
library(shinydashboard)
library(dashboardthemes)
library(janitor)



# Functions ----------------------------------------------------------------------

volcanoplot <- function(res, lfcthresh = 2, sigthresh = 0.05, main = "Volcano Plot", legendpos = "bottomright", labelsig = TRUE, textcx = 1, ...) {
    with(res, plot(log2FoldChange, -log10(pvalue), pch = 20, main = main, ...))
    with(subset(res, padj < sigthresh), points(log2FoldChange, -log10(pvalue), pch = 20, col = "red", ...))
    with(subset(res, abs(log2FoldChange) > lfcthresh), points(log2FoldChange, -log10(pvalue), pch = 20, col = "orange", ...))
    with(subset(res, padj < sigthresh & abs(log2FoldChange) > lfcthresh), points(log2FoldChange, -log10(pvalue), pch = 20, col = "green", ...))
    legend(legendpos, xjust = 1, yjust = 1, legend = c(paste("FDR<", sigthresh, sep = ""), paste("|LogFC|>", lfcthresh, sep = ""), "both"), pch = 20, col = c("red", "orange", "green"))
}


# Preliminary code ---------------------------------------------------------------
# update max upload size to 50mB
options(shiny.maxRequestSize = 50*1024^2)
enr_sources <- c("GO", "KEGG", "REAC", "TF", "MIRNA", "CORUM", "HP", "HPA", "WP")

# This loads a DESeq object prealably saved (needs to be saved)
# and avoids calculations
testing <- FALSE 

# This permits toggle between tximeta (TRUE) and tximport (FALSE) (does not work)
txi_met_chosen <- TRUE
