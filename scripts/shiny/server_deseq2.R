# Validate the input of sft_thres
iv$add_rule("sft_thres", sv_between(left = 1, right = 30))
iv$enable()

# enforces that base condition and compare condition are different
observe(iv$add_rule("compare_cond",
            sv_not_equal(rhs = input$base_cond,
                                         message_fmt = "This condition should be different from the base one")))



observe({
  req(my_values$coldata)
  updateSelectInput(session = session,
                    inputId = "variables",
                    choices = colnames(my_values$coldata))
})

observe({
  req(my_values$coldata)
  updateSelectInput(session = session,
                    inputId = "deseq_var",
                    choices = colnames(my_values$coldata))
})


observe({
  req(my_values$coldata)
  req(input$deseq_var)
  my_values$variable_chosen <- my_values$coldata %>% dplyr::pull(input$deseq_var) %>% as.factor()
})

observe({
  req(my_values$coldata)
  req(input$deseq_var)
  updateSelectInput(session = session,
                    inputId = "base_cond",
                    choices = unique(my_values$coldata[,input$deseq_var]))
})

observe({
  req(my_values$coldata)
  req(input$deseq_var)
  updateSelectInput(session = session,
                    inputId = "compare_cond",
                    choices = unique(my_values$coldata[,input$deseq_var]))
})

#The table of genes, displayed on DEG panel
gene_table <- reactive({
  my_values$res[order(my_values$res$padj), ] %>% 
    as.data.frame() %>% 
    filter(padj < input$pval_cutoff, 
           log2FoldChange > input$lfc_cutoff |
             log2FoldChange < -input$lfc_cutoff) %>%
    #Significant digits
    mutate(dplyr::across(where(is.numeric), signif, 3))
})

#The button to run DESeqDataSet, DESeq, results
dds <- eventReactive(input$execute_d, {
  # relevel based on the user selection
  my_values$variable_chosen %<>% relevel(input$base_cond)
  # my_values$coldata[, input$deseq_var] %<>% relevel(input$base_cond) #usefulness?
  withProgress(message = "Running DESeq2", {
    if(input$snakemake) {
      withProgress(message = "Loading data...", {
        if(!txi_met_chosen) {
          tmp <- DESeqDataSetFromTximport(txi,
                                          my_values$coldata,
                                          paste("~ ",
                                                paste(req(input$variables),
                                                      collapse = " + ")) %>%
                                            as.formula())
        } else {
          tmp <- DESeqDataSet(gse(), 
                              design = paste("~ ",
                                             paste(req(input$variables),
                                                   collapse = " + ")) %>% 
                                as.formula())
        }
      })
    }
    else {
      if(testing) {
        load(file = "./dds.RData")
        return(tmp2)
      }
      else {
      tmp <- DESeqDataSetFromMatrix(req(my_values$counts),
                                    my_values$coldata,
                                    design = paste("~ ",
                                                   paste(req(input$variables),
                                                         collapse = " + ")) %>% 
                                      as.formula())
      }
    }
    #filtering
    keep <- rowSums(counts(tmp)) > 1
    tmp2 <- DESeq(tmp[keep,])
    # save(tmp2, file = "./dds.RData")
    return(tmp2)
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
  rownames(sampleDistMatrix) <- my_values$coldata[, "names"]
  colnames(sampleDistMatrix) <- my_values$coldata[, "names"]
  
  colors <- colorRampPalette(rev(brewer.pal(9, "Purples")) )(255)
  pheatmap(sampleDistMatrix,
           show_rownames = FALSE,
           show_colnames = TRUE,
           annotation_col = my_values$coldata %>%
             select("names", input$variables) %>%
             tibble::column_to_rownames("names"),
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
}, res = 96, height = 600)

output$pca <- renderPlot({
  req(rld_vst())
  plotPCA(rld_vst(), intgroup = input$variables, ntop = 1000)
}, res = 96, height = 600)

output$ma <- renderPlot({
  req(my_values$res)
  plotMA(my_values$res, alpha = 0.05)
}, res = 96)


output$volcano <- renderPlot({
  req(my_values$res)
  volcanoplot(my_values$res,
              lfcthresh=1,
              sigthresh=0.05,
              textcx=.8,
              xlim=c(-2.3, 2))
}, res = 96, height = 600)

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
  
  return(req(plot_counts()))
}, res = 96)

plot_counts <- reactive({
  req(input$sel_gene)
  req(my_values$res)
   
  d <- plotCounts(dds(),
                  gene=which(my_values$res$symbol == input$sel_gene),
                  intgroup = input$deseq_var,
                  returnData=TRUE)
  
  plot_fin <- ggplot(data = d %>% filter(my_values$variable_chosen %in% input$condition_plot),
                     aes_string(x = input$deseq_var, y = "count")) + 
    geom_point(position = position_jitter(w=0.1,h=0)) + 
    ggtitle(paste("Plot of ", input$sel_gene, " counts")) + 
    xlab("Modality") +
    ylab("Counts") +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12))
  if(input$plot == "base") {
    plot_fin <- plot_fin + scale_y_log10()
  }
  if(input$plot == "barplot"){
    plot_fin <- plot_fin + stat_summary(aes_string(fill = input$deseq_var),
                                            fun = mean,
                                            geom = "bar",
                                            alpha = input$alpha) +
      stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.05)
  }
  if(input$plot == "boxplot") {
   plot_fin <- plot_fin + geom_boxplot(aes_string(fill = input$deseq_var))
  }
  plot_fin
})


output$down_plot <- downloadHandler(
  filename = function() {
    paste(req(input$sel_gene), ".svg", sep = "")
  },
  content = function(file) {
    ggsave(plot = req(plot_counts()), filename = file, device = "svg")
  }
)