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
library(dplyr)
library(janitor)
library(scales)



# Preliminary code ---------------------------------------------------------------
# update max upload size to 50mB
options(shiny.maxRequestSize = 50*1024^2)
enr_sources <- c("GO", "KEGG", "REAC", "TF", "MIRNA", "CORUM", "HP", "HPA", "WP")

# This loads a DESeq object prealably saved (needs to be saved)
# and avoids calculations
testing <- FALSE 

# This permits toggle between tximeta (TRUE) and tximport (FALSE does not work)
txi_met_chosen <- TRUE
