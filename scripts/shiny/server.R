server <- function(input, output, session) {
  
  #An object to store everything
  my_values <- reactiveValues()
  
  #The table of genes, displayed on DEG panel
  gene_table <- eventReactive(my_values$res, {
    my_values$res[order(my_values$res$padj), ] %>% 
      as.data.frame() %>% 
      filter(padj < input$pval_cutoff, 
             log2FoldChange > input$lfc_cutoff |
             log2FoldChange < -input$lfc_cutoff) %>%
      #Significant digits
      mutate(across(where(is.numeric), signif, 3))
  })

  #The button to run DESeqDataSet, DESeq, results
  observeEvent(input$execute_d, {
    req(input$variables)
    req(input$compare_cond)
    gse$condition %<>% relevel(input$base_cond)
    my_values$dds <- DESeqDataSet(gse, 
                                  design = paste("~ ",
                                                 paste(input$variables,
                                                       collapse = " + ")) %>% 
                                    as.formula())
    #filtering
    keep <- rowSums(counts(my_values$dds)) > 1 
    my_values$dds <- DESeq(my_values$dds[keep,])
    
    
    my_values$rld <- rlogTransformation(my_values$dds)
    
    #temporary first element of contrast /!\
    my_values$res <- results(my_values$dds,
                             contrast = c(input$variables[1],
                                          input$compare_cond, input$base_cond))
    
    #Adding gene names
    # /!\/!\ should look at duplicates and is.na to look for lack of info
    ens.str <- substr(rownames(my_values$res), 1, 15)
    my_values$res$symbol <- mapIds(org.Hs.eg.db,
                                   keys=ens.str,
                                   column="SYMBOL",
                                   keytype="ENSEMBL",
                                   multiVals="first")
    
  })
  
  observeEvent(input$explore_w, {
    #using bias corrected counts (without an offset)
    raw_counts <- assay(gse_scaled, "counts") %>% round()
    # Remove non-expressed genes
    raw_counts <- raw_counts[rowSums(raw_counts) > 0,]
    
    # Calculate upper quartile value (not taking into account 0s by sample)
    quantile_expressed <- apply(raw_counts, 
                                2, 
                                function(x){quantile(x[x>0], 0.75)})
    # Divide each column by its upper quartile value
    my_values$counts_norm <- t(raw_counts)/quantile_expressed
  })
  
  observeEvent(input$rm_sample, {
    my_values$counts_norm <- 
        my_values$counts_norm[!(rownames(my_values$counts_norm) %in%
                                input$rm_sample),]
  })
  
  output$dist <- renderPlot({
    req(my_values$rld)
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
    req(my_values$rld)
    plotPCA(my_values$rld)
  })
  
  output$ma <- renderPlot({
    req(my_values$res)
    plotMA(my_values$res, alpha = 0.05)
  })
  
  
  output$volcano <- renderPlot({
    req(my_values$res)
    volcanoplot(my_values$res,
                lfcthresh=1,
                sigthresh=0.05,
                textcx=.8,
                xlim=c(-2.3, 2))
  })
  
  output$nb_genes <- renderText({
    req(my_values$res)
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
    
    paste("There are ", up,
          " significantly upregulated genes", " and there are ",
          down, "significantly downregulated genes", "at a p-value of", 
          input$pval_cutoff, " and a LFC of ", input$lfc_cutoff)
  })
  
  output$genes <- DT::renderDataTable({
    req(my_values$res)
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
      scale_y_log10()
  })
  
  observeEvent(my_values$counts_norm, {
    # Update the maximum of gene one can keep
    # default value of 8 000 genes not to crush computer resources
    m <- round(8000/ncol(my_values$counts_norm), 1)
    updateSliderInput(session,
                      "percent_g",
                      max = m)
  })
  
  observeEvent(input$percent_g, {
    req(my_values$counts_norm)
    my_values$counts_filt <- filter_low_var(my_values$counts_norm,
                                            pct = input$percent_g,
                                            type = "median")
  })
  
  output$outliers <- renderPlot({
    req(input$explore_w)
    sampleTree <- hclust(dist(my_values$counts_norm), method = "average")
    plot(sampleTree, main = "Sample clustering to detect outliers",
         sub="", xlab="", cex.lab = 1.5,
         cex.axis = 1.5, cex.main = 2)
  })
  
  output$warning <- renderText({
    req(my_values$counts_norm)
    if(nrow(my_values$counts_norm) < 20) {
      paste("\n", "Warning : The number of samples is too low for WGCNA analyis.", "\n")
    }
  })
  
  output$threshold <- renderText({
    req(my_values$counts_filt)
    
    # Rerun when button is clicked
    # But do not run when not clicked
    if(input$update_sft == 0)
      return()
    
    pwr_vec <- c(1:9, seq(10, 30, by = 2)) 
    withProgress(message = "Calculating", {
      isolate(sft <- WGCNA::pickSoftThreshold(my_values$counts_filt,
                                      powerVector = pwr_vec,
                                      corFnc = WGCNA::cor,
                                      corOptions = list(method = "spearman"),
                                      networkType = "signed hybrid"))
    })
    if(is.na(sft$powerEstimate)) {
      max_pow <- max(sft$fitIndices[, "Power"])
      paste("WGCNA's pickSoftThreshold did not find an optimal power value.",
            "Max SFT.R.sq :  ",
            round(sft$fitIndices[sft$fitIndices$Power == max_pow, "SFT.R.sq"],2),
            " Power : ", max_pow) 
    } else {
      paste("WGCNA's pickSoftThreshold recommends a power threshold of ",
            sft$powerEstimate,
            "which corresponds to a R squared for SFT of ",
            round(sft$fitIndices[sft$fitIndices$Power == sft$powerEstimate, "SFT.R.sq"],2)) 
    }
  })
  
  
}
