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

ui <- fluidPage(
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
    
    tabPanel("Most DEGs",
             sidebarLayout(
               
               sidebarPanel(
                 sliderInput(inputId = "n_genes",
                             label = "Number of most differentially expressed genes shown : ",
                             min = 1,
                             max = 1000,
                             value = 50,
                             step = 10,
                             round = 1)
               ),
               
               mainPanel(
                 DT::dataTableOutput(outputId = "genes")
               )
             )
    )
    
    # tabPanel("Volcano plot", plotOutput(outputId = "volcano"))
  )
)


# Server -------------------------------------------------------------------------
server <- function(input, output, session) {
  
  my_values <- reactiveValues()
  
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
    my_values$res <- results(my_values$dds,
                             contrast = c(input$variables[1], input$compare_cond, input$base_cond))
    
    
  })
  
  
  output$pca <- renderPlot({
    plotPCA(my_values$rld)
  })
  
  output$genes <- DT::renderDataTable({
    head(my_values$res[order(my_values$res$padj), ], input$n_genes) %>% 
      as.data.frame() %>% 
      select(baseMean, log2FoldChange, padj) %>%
      #Significant digits
      mutate(across(everything(), signif, 3))
  })
  
  # 
  # output$volcano <- renderPlot(
  #   volcanoplot(res(), lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
  # )
}

# Shiny Object -------------------------------------------------------------------
shinyApp(ui = ui, server = server)