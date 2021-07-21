server <- function(input, output, session) {
  
  #An object to store everything
  my_values <- reactiveValues()
  
  #The table of genes, displayed on DEG panel
  gene_table <- reactive({
    my_values$res[order(my_values$res$padj), ] %>% 
      as.data.frame() %>% #or assay()?
      filter(padj < input$pval_cutoff, 
             log2FoldChange > input$lfc_cutoff | log2FoldChange < -input$lfc_cutoff) %>%
      #Significant digits
      mutate(across(where(is.numeric), signif, 3))
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
    
    #Adding gene names
    ens.str <- substr(rownames(my_values$res), 1, 15)
    my_values$res$symbol <- mapIds(org.Hs.eg.db,
                                   keys=ens.str,
                                   column="SYMBOL",
                                   keytype="ENSEMBL",
                                   multiVals="first")
    
  })
  
  
  output$dist <- renderPlot({
    sampleDists <- my_values$rld %>%
      assay() %>%
      t() %>%
      dist()
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- my_values$rld$condition
    colnames(sampleDistMatrix) <- NULL
    
    colors <- colorRampPalette( rev(brewer.pal(9, "Purples")) )(255)
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             col=colors)
  })
  
  output$pca <- renderPlot({
    plotPCA(my_values$rld)
  })
  
  output$ma <- renderPlot({
    plotMA(my_values$res, alpha = 0.05)
  })
  
  
  output$volcano <- renderPlot({
    volcanoplot(my_values$res, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
  })
  
  output$nb_genes <- renderText({
    req(input$pval_cutoff)
    req(input$lfc_cutoff)
    
    up <- my_values$res %>% 
      as.data.frame() %>% #Or assay()
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
  
  observe({
    updateSelectizeInput(session,
                         inputId = "sel_gene", 
                         label = "Select which genes to plot :",
                         choices = as.vector(my_values$res$symbol), 
                         server = TRUE)
  })
  
  
  output$plot_gene <- renderPlot({
    req(input$sel_gene)
    d <- plotCounts(my_values$dds,
                    gene=which(my_values$res$symbol == input$sel_gene),
                    returnData=TRUE)
    ggplot(d, aes(x=condition, y=count)) + 
      geom_point(position=position_jitter(w=0.1,h=0)) + 
      scale_y_log10(breaks=c(25,100,400))
  })
  
  
  
}