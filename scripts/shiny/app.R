# load libraries -----------------------------------------------------------------
library(shiny)
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
library(shinythemes)

# Functions ----------------------------------------------------------------------
# volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
#   with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
#   with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
#   with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
#   with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
#   legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
# }


# Preliminary code ---------------------------------------------------------------
load("txidata.RData")



# UI -----------------------------------------------------------------------------

ui <- fluidPage(theme = shinytheme("cosmo"),
  navbarPage(
    title = "My little app",
    tabPanel("Setup",
             #Choose the experimental design
             selectInput(inputId = "variables",
                         label = "Choose the variables of the experimental design :",
                         choices = colnames(coldata),
                         multiple = T),
             
             # Choose base level for condition
             selectInput(inputId = "base_cond",
                         label = "Choose the base condition :",
                         choices = levels(gse$condition)),
             
             # Choose base level for condition
             selectInput(inputId = "compare_cond",
                         label = "Choose the condition to compare with:",
                         choices = levels(gse$condition)),
             
             actionButton(inputId = "execute",
                          label = "Run DESeq2")
    ),
    tabPanel("PCA", plotOutput(outputId = "pca")),
    
    tabPanel("MAplot", plotOutput(outputId = "ma")),
    
    tabPanel("Most DEGs",
             sidebarLayout(
               
               sidebarPanel(
                 numericInput(inputId = "pval_cutoff",
                             label = "Enter the maximum p-value :",
                             value = 0.05,
                             min = 0,
                             max = 1,
                             step = .05),
                 
                 numericInput(inputId = "lfc_cutoff",
                              label = "Enter the minimum (absolute) logFoldChange :",
                              value = 1,
                              min = 0)
               ),
               
               mainPanel(
                 
                 textOutput(outputId = "nb_genes"),
                 br(),
                 DT::dataTableOutput(outputId = "genes"),

		 downloadButton(outputId = "download",
                                label = "Download table")
               )
             )
    )
    
    # tabPanel("Volcano plot", plotOutput(outputId = "volcano"))
  )
)


# Server -------------------------------------------------------------------------
server <- function(input, output, session) {
  
  #An object to store everything
  my_values <- reactiveValues()
  
  #The table of genes, displayed on DEG panel
  gene_table <- reactive({
    my_values$res[order(my_values$res$padj), ] %>% 
      as.data.frame() %>% 
      select(baseMean, log2FoldChange, padj) %>%
      filter(padj < input$pval_cutoff, log2FoldChange > input$lfc_cutoff | log2FoldChange < -input$lfc_cutoff) %>%
      #Significant digits
      mutate(across(everything(), signif, 3))
  })
  
  #The button to run DESeqDataSet, DESeq, results
  observeEvent(input$execute, {
    req(input$variables)
    req(input$compare_cond)
    gse$condition %<>% relevel(input$base_cond)
    my_values$dds <- DESeqDataSet(gse, 
                                  design = paste("~ ", paste(input$variables, collapse = " + ")) %>% 
                                    as.formula())
    #filtering
    keep <- rowSums(counts(my_values$dds)) > 1 
    my_values$dds <- DESeq(my_values$dds[keep,])
    
    
    my_values$rld <- rlogTransformation(my_values$dds)
    
    #temporary first element of contrast /!\
    my_values$res <- results(my_values$dds,
                             contrast = c(input$variables[1], input$compare_cond, input$base_cond))
    
    
  })
  
  
  output$pca <- renderPlot({
    plotPCA(my_values$rld)
  })
  
  output$ma <- renderPlot({
    plotMA(my_values$res, alpha = 0.05)
  })
  
  output$nb_genes <- renderText({
    req(input$pval_cutoff)
    req(input$lfc_cutoff)
    
    up <- my_values$res %>% 
      as.data.frame() %>%
      filter(padj < input$pval_cutoff, log2FoldChange > input$lfc_cutoff) %>%
      nrow()
    
    down <- my_values$res %>% 
      as.data.frame() %>%
      filter(padj < input$pval_cutoff, log2FoldChange < -input$lfc_cutoff) %>%
      nrow()
    
    paste("There are ", up, " significantly upregulated genes", " and there are ",
          down, "significantly downregulated genes", "at a p-value of", 
          input$pval_cutoff, " and a LFC of ", input$lfc_cutoff)
  })
  
  output$genes <- DT::renderDataTable({
    req(input$pval_cutoff)
    req(input$lfc_cutoff)
    gene_table()
  })
  
  output$download <- downloadHandler(
    filename = function() {
      paste("genes", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(gene_table(), file)
    }
  )
  
  # 
  # output$volcano <- renderPlot(
  #   volcanoplot(res(), lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
  # )
}

# Shiny Object -------------------------------------------------------------------
shinyApp(ui = ui, server = server)