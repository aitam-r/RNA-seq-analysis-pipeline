# Validate the input of sft_thres
iv$add_rule("sft_thres", sv_between(left = 1, right = 30))
iv$enable()

# Object to store reactives that can change
my_values <- reactiveValues(res = NULL, counts_norm = NULL, counts_filt = NULL)


#The table of genes, displayed on DEG panel
gene_table <- reactive({
  my_values$res()[order(my_values$res()$padj), ] %>% 
    as.data.frame() %>% 
    filter(padj < input$pval_cutoff, 
           log2FoldChange > input$lfc_cutoff |
             log2FoldChange < -input$lfc_cutoff) %>%
    #Significant digits
    mutate(across(where(is.numeric), signif, 3))
})

#The button to run DESeqDataSet, DESeq, results
dds <- eventReactive(input$execute_d, {
  req(input$variables)
  req(input$compare_cond)
  gse$condition %<>% relevel(input$base_cond)
  withProgress(message = "Running DESeq2", {
    tmp <- DESeqDataSet(gse, 
                        design = paste("~ ",
                                       paste(input$variables,
                                             collapse = " + ")) %>% 
                          as.formula())
    #filtering
    keep <- rowSums(counts(tmp)) > 1 
    DESeq(tmp[keep,])
  })
})

rld <- eventReactive(dds(), {
  rlog(dds())
})

my_values$res <- eventReactive(dds(), {
  #temporary first element of contrast /!\
  tmp <- results(dds(),
                 contrast = c(input$variables[1],
                              input$compare_cond, input$base_cond))
  
  #Adding gene names
  # /!\/!\ should look at duplicates and is.na to look for lack of info
  
  ens.str <- substr(rownames(tmp), 1, 15)
  tmp$symbol <- mapIds(org.Hs.eg.db,
                       keys=ens.str,
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
  tmp
})



output$dist <- renderPlot({
  req(rld())
  sampleDists <- rld() %>%
    assay() %>%
    t() %>%
    dist()
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- rld()$condition
  colnames(sampleDistMatrix) <- NULL
  
  colors <- colorRampPalette(rev(brewer.pal(9, "Purples")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
})

output$pca <- renderPlot({
  req(rld())
  plotPCA(rld())
})

output$ma <- renderPlot({
  req(my_values$res())
  plotMA(my_values$res(), alpha = 0.05)
})


output$volcano <- renderPlot({
  req(my_values$res())
  volcanoplot(my_values$res(),
              lfcthresh=1,
              sigthresh=0.05,
              textcx=.8,
              xlim=c(-2.3, 2))
})

# outputs the number of genes without names
output$genes_na <- renderUI({
  req(gene_table())
  nb_na <- length(which(is.na(gene_table()[,"symbol"])))
  HTML(paste0("<p> There are ", nb_na, " genes without an affected gene name"))
})

output$nb_genes <- renderText({
  req(my_values$res())
  req(input$pval_cutoff)
  req(input$lfc_cutoff)
  
  # upregulated genes
  up <- my_values$res() %>% 
    as.data.frame() %>% 
    filter(padj < input$pval_cutoff, log2FoldChange > input$lfc_cutoff) %>%
    nrow()
  
  # downregulated genes
  down <- my_values$res() %>% 
    as.data.frame() %>%
    filter(padj < input$pval_cutoff, log2FoldChange < -input$lfc_cutoff) %>%
    nrow()
  
  paste("There are ", up,
        " significantly upregulated genes", " and there are ",
        down, "significantly downregulated genes", "at a p-value of", 
        input$pval_cutoff, " and a LFC of ", input$lfc_cutoff)
})

output$genes <- DT::renderDataTable({
  req(gene_table())
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

observe({
  updateSelectizeInput(session,
                       inputId = "sel_gene", 
                       label = "Select which genes to plot :",
                       choices = as.vector(my_values$res()$symbol), 
                       server = TRUE)
})

observe({
  updateCheckboxGroupInput(session,
                           inputId = "condition_plot",
                           selected = c(input$base_cond, input$compare_cond))
})

output$plot_gene <- renderPlot({
  req(input$sel_gene)
  d <- plotCounts(dds(),
                  gene=which(my_values$res()$symbol == input$sel_gene),
                  returnData=TRUE)
  if(input$barplot) {
  ggplot(data = d %>% filter(condition %in% input$condition_plot),
              aes(x=condition, y=count)) + 
    geom_point(position=position_jitter(w=0.1,h=0)) + 
    scale_y_log10() +
    # Many assumptions here, but at least not normality (mean_cl_boot)
    stat_summary(fun = mean, geom = "bar", alpha = 0.2) +
       stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.05)
  } else {
    ggplot(data = d %>% filter(condition %in% input$condition_plot),
              aes(x=condition, y=count)) + 
    geom_point(position=position_jitter(w=0.1,h=0)) + 
    scale_y_log10()
  }
})
