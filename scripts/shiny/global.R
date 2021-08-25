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



# Functions ----------------------------------------------------------------------

volcanoplot <- function(res, lfcthresh = 2, sigthresh = 0.05, main = "Volcano Plot", legendpos = "bottomright", labelsig = TRUE, textcx = 1, ...) {
    with(res, plot(log2FoldChange, -log10(pvalue), pch = 20, main = main, ...))
    with(subset(res, padj < sigthresh), points(log2FoldChange, -log10(pvalue), pch = 20, col = "red", ...))
    with(subset(res, abs(log2FoldChange) > lfcthresh), points(log2FoldChange, -log10(pvalue), pch = 20, col = "orange", ...))
    with(subset(res, padj < sigthresh & abs(log2FoldChange) > lfcthresh), points(log2FoldChange, -log10(pvalue), pch = 20, col = "green", ...))
    legend(legendpos, xjust = 1, yjust = 1, legend = c(paste("FDR<", sigthresh, sep = ""), paste("|LogFC|>", lfcthresh, sep = ""), "both"), pch = 20, col = c("red", "orange", "green"))
}
get_sub_clusters <- function(network, seq_k = seq_len(15), fit_plot = TRUE,
                             ...) {
  # Checking args
  is_network(network)
  if(!is.numeric(seq_k))
    stop("seq_k must be a numeric vector")
  if(length(seq_k) < 2)
    stop("seq_k must at leat have two power values to test")
  if (any(lapply(seq_k, function(x) x < 1 | x %% 1 != 0) %>% unlist))
    stop("seq_k must contain only whole numbers >= 1")
  if(!is.logical(fit_plot))
    stop("fit_plot must be a boolean")
  # Computing the k-medoid
  list_k_tests <- lapply(seq_k, function(k) {
    k_res <- cluster::pam(network, k, diss = TRUE, ...)
    return(k_res)
  }) %>% stats::setNames(paste0("k_", seq_k))
  
  # Summarizing needed results for plot and returning optimal k
  df_k_tests <- list_k_tests %>%
    lapply(`[[`, "silinfo") %>%
    lapply(`[[`, "avg.width") %>%
    unlist %>%
    data.frame(k = seq_k[-1], avg_sil_width = .)
  
  # If asked, plotting the silhouette coefficient against k tested
  if (fit_plot) {
    graphics::plot(df_k_tests, type="b", pch = 19,
                   frame = FALSE, xlab="Number of clusters K",
                   ylab="Average silhouette width")
  }
  
  # Isolating optimal k regarding the silhouette coefficient
  opti_k <- df_k_tests %>%
    dplyr::filter(avg_sil_width == max(avg_sil_width)) %>%
    dplyr::select(k) %>% as.numeric()
  
  # Formatting the table to return the gene id association to sub_module
  clusters <- list_k_tests[[opti_k]]$clustering %>%
    data.frame(sub_module = .) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(sub_module = as.character(sub_module))
  
  return(clusters)
}


# Preliminary code ---------------------------------------------------------------
# update max upload size to 50mB
options(shiny.maxRequestSize = 50*1024^2)
enr_sources <- c("GO", "KEGG", "REAC", "TF", "MIRNA", "CORUM", "HP", "HPA", "WP")
testing <- TRUE 
txi_met_chosen <- TRUE
