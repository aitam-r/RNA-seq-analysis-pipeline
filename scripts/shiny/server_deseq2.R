# Validate the input of sft_thres
iv$add_rule("sft_thres", sv_between(left = 1, right = 30))
iv$enable()

# enforces that base condition and compare condition are different
observe(iv$add_rule("compare_cond",
            sv_not_equal(rhs = input$base_cond,
                                         message_fmt = "This condition should be different from the base one")))



observe({
  req(coldata())
  updateSelectInput(session = session,
                    inputId = "variables",
                    choices = colnames(coldata()))
})

observe({
  req(coldata())
  updateSelectInput(session = session,
                    inputId = "deseq_var",
                    choices = colnames(coldata()))
})


observe({
  req(coldata())
  req(input$deseq_var)
  my_values$variable_chosen <- coldata() %>% pull(input$deseq_var) %>% as.factor()
})

observe({
  req(coldata())
  req(input$deseq_var)
  updateSelectInput(session = session,
                    inputId = "base_cond",
                    choices = unique(coldata()[,input$deseq_var]))
})

observe({
  req(coldata())
  req(input$deseq_var)
  updateSelectInput(session = session,
                    inputId = "compare_cond",
                    choices = unique(coldata()[,input$deseq_var]))
})

#The table of genes, displayed on DEG panel
gene_table <- reactive({
  my_values$res[order(my_values$res$padj), ] %>% 
    as.data.frame() %>% 
    filter(padj < input$pval_cutoff, 
           log2FoldChange > input$lfc_cutoff |
             log2FoldChange < -input$lfc_cutoff) %>%
    #Significant digits
    mutate(across(where(is.numeric), signif, 3))
})

#The button to run DESeqDataSet, DESeq, results
dds <- eventReactive(input$execute_d, {
  
  # relevel based on the user selection
  my_values$variable_chosen %<>% relevel(input$base_cond)
  
  withProgress(message = "Running DESeq2", {
    if(input$snakemake) {
      withProgress(message = "Loading data...", {
        #fonction_load()
        load("~/Documents/TAAAAF/Stage/liver/scripts/shiny/txidata.RData")
      })
      tmp <- DESeqDataSet(gse, 
                          design = paste("~ ",
                                         paste(req(input$variables),
                                               collapse = " + ")) %>% 
                            as.formula())
    }
    else {
     
      
      tmp <- DESeqDataSetFromMatrix(req(my_values$counts),
                                    coldata(),
                                    design = paste("~ ",
                                                   paste(req(input$variables),
                                                         collapse = " + ")) %>% 
                                      as.formula())
    }
    #filtering
    keep <- rowSums(counts(tmp)) > 1 
    DESeq(tmp[keep,])
  })
})

rld_vst <- eventReactive(dds(), {
  req(input$r_v)
  if(input$r_v == "rlog"){
  withProgress(message = "Calculating rlog", {
  rlog(dds())
  })
  } else {
    withProgress(message = "Calculating vst", {
      vst(dds())
    })
  }
})

observeEvent(dds(), {
  withProgress(message = "Getting results", { 
  my_values$res <- results(dds(),
                 contrast = c(input$deseq_var,
                              input$compare_cond, input$base_cond))
  })
  
  #Adding gene names
  # /!\/!\ should look at duplicates and is.na to look for lack of info
  
  ens.str <- substr(rownames(my_values$res), 1, 15)
  my_values$res$symbol <- mapIds(org.Hs.eg.db,
                       keys=ens.str,
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
})



output$dist <- renderPlot({
  req(rld_vst())
  sampleDists <- rld_vst() %>%
    assay() %>%
    t() %>%
    dist()
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- rld_vst()$condition
  colnames(sampleDistMatrix) <- NULL
  
  colors <- colorRampPalette(rev(brewer.pal(9, "Purples")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
})

output$pca <- renderPlot({
  req(rld_vst())
  plotPCA(rld_vst(), intgroup = input$variables)
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

# outputs the number of genes without names
output$genes_na <- renderUI({
  req(gene_table())
  nb_na <- length(which(is.na(gene_table()[,"symbol"])))
  HTML(paste0("<p> There are ", nb_na, " genes without an affected gene name"))
})

output$nb_genes <- renderText({
  req(my_values$res)
  req(input$pval_cutoff)
  req(input$lfc_cutoff)
  
  # upregulated genes
  up <- my_values$res %>% 
    as.data.frame() %>% 
    filter(padj < input$pval_cutoff, log2FoldChange > input$lfc_cutoff) %>%
    nrow()
  
  # downregulated genes
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

observeEvent(my_values$res, {
  updateSelectizeInput(session,
                       inputId = "sel_gene", 
                       label = "Select which genes to plot :",
                       choices = as.vector(my_values$res$symbol),
                       server = TRUE)
})

observe({
  updateCheckboxGroupInput(session,
                           inputId = "condition_plot",
                           choices = levels(my_values$variable_chosen),
                           selected = c(input$base_cond, input$compare_cond))
})

output$plot_gene <- renderPlot({
  req(input$sel_gene)
  req(my_values$res)
  
  # find the minimum of replicates by condition (for barplot transparency)
  min_replicates <- count(my_values$variable_chosen %>% as.data.frame(), 
                          my_values$variable_chosen)[,"n"] %>%
    min()
  
  d <- plotCounts(dds(),
                  gene=which(my_values$res$symbol == input$sel_gene),
                  intgroup = input$deseq_var,
                  returnData=TRUE)
  
  levels_variable <- levels(my_values$variable_chosen)
  
  if(input$barplot) {
  ggplot(data = d %>% filter(my_values$variable_chosen %in% input$condition_plot),
              aes(x = my_values$variable_chosen, y = count)) + 
    geom_point(position = position_jitter(w=0.1,h=0)) + 
    scale_y_log10() +
    # Many assumptions here, but at least not normality (mean_cl_boot)
    stat_summary(fun = mean, geom = "bar", alpha = 0.05*min_replicates) +
       stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.05)
  } else {
    ggplot(data = d %>% filter(my_values$variable_chosen %in% input$condition_plot),
              aes(x = my_values$variable_chosen, y = count)) + 
    geom_point(position = position_jitter(w=0.1,h=0)) + 
    scale_y_log10()
  }
})
