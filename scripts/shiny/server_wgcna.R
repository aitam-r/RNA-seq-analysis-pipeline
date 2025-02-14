# Calculs ------

# Obtains normalized counts, depending on the provenance
observe({
  req(input$upload_samp_tab)
  
  # If the counts are from the pipeline
  if(input$snakemake) {
  req(my_values$se) 
    if(!txi_met_chosen){ # From tximport
      txi_scaled <- tximport(req(my_values$coldata$files),
                             type = "salmon",
                             tx2gene = tx2gene(),
                             countsFromAbundance = "lengthScaledTPM")
      raw_counts <- txi_scaled$counts %>% round()
    } else{ # From tximeta
      gse_scaled <- summarizeToGene(my_values$se, countsFromAbundance = "lengthScaledTPM")
      raw_counts <- assay(gse_scaled, "counts") %>% round()
    }
  } else { # If counts are given in the shiny app
    req(my_values$counts)
    raw_counts <- isolate(my_values$counts) %>% round()
  }
  # Remove non-expressed genes
  raw_counts <- raw_counts[rowSums(raw_counts) > 0,]
  
  # Calculate upper quartile value (not taking into account 0s by sample)
  quantile_expressed <- apply(raw_counts, 
                              2, 
                              function(x){quantile(x[x>0], 0.75)})
  # Divide each column by its upper quartile value
  my_values$counts_norm <- t(raw_counts)/quantile_expressed
})

# Check the minimum number of samples for WGCNA analysis
output$warning <- renderText({
  req(my_values$counts_norm)
  if(nrow(my_values$counts_norm) < 20) {
    paste("\n", "Warning : The number of samples is too low for WGCNA analyis.", "\n")
  }
})

# Classic hclust plot to screen for outliers
output$outliers <- renderPlot({
  req(input$explore_w)
  sampleTree <- hclust(dist(my_values$counts_norm), method = "average")
  plot(sampleTree, main = "Sample clustering to detect outliers",
       sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
})

# Remove a sample
observeEvent(input$rm_sample, {
  my_values$counts_norm <- 
    my_values$counts_norm[!(rownames(my_values$counts_norm) %in%
                              input$rm_sample),]
})

observe({
  req(my_values$coldata)
  updateSelectizeInput(inputId = "rm_sample",
                       choices = my_values$coldata[, "names"],
                       options = list(maxItems = nrow(my_values$coldata) - 3))
})


net <- eventReactive(input$build, {
  withProgress(message = "Building network", {
    build_net(my_values$counts_filt,
              cor_func = "spearman",
              power_value = input$sft_thres,
              network_type = isolate(input$type_net),
              n_threads = 1)
  })
})


modules <- eventReactive(net(), {
  withProgress(message = "Detecting Modules", {
    detect_modules(my_values$counts_filt,
                   net()$network,
                   detailled_result = T)
  })
})

modules_enriched <- eventReactive({
  modules()
  input$enrich}, { 
    withProgress(message = "Performing Enrichment", {
      bio_enrich(modules()$modules,
                 organism = "hsapiens",
                 sources = input$sources)
    })
  })

observeEvent(my_values$counts_norm, {
  # Update the maximum number of gene one can keep
  m <- req(ncol(my_values$counts_norm))
  updateSliderInput(session,
                    "number_g",
                    max = m)
})

observe({
  req(my_values$counts_norm)
  req(input$number_g)
  my_values$counts_filt <- filter_low_var(my_values$counts_norm,
                 pct = input$number_g/ncol(my_values$counts_norm),
                 type = "median")
})

# tom_diss <- eventReactive(net(), {
#   WGCNA::TOMdist(net()$network)
# })
# 
# 
# output$tomplot <- renderPlot({
#   WGCNA::TOMplot(tom_diss(), fastcluster::hclust(as.dist(tom_diss()), method = "average"))
# })


sft <- eventReactive({input$update_sft}, {
  withProgress(message = "Calculating", {
    pwr_vec <- c(1:9, seq(10, 30, by = 2))
    WGCNA::pickSoftThreshold(req(my_values$counts_filt),
                             powerVector = pwr_vec,
                             corFnc = WGCNA::cor,
                             corOptions = list(method = "spearman"),
                             networkType = input$type_net)
  })
})

output$threshold <- renderText({
  req(sft())
  if(is.na(sft()$powerEstimate)) {
    max_pow <- max(sft()$fitIndices[, "Power"])
    paste("WGCNA's pickSoftThreshold did not find an optimal power value.",
          "Max SFT.R.sq :  ",
          round(sft()$fitIndices[sft()$fitIndices$Power == max_pow, "SFT.R.sq"],2),
          " Power : ", max_pow) 
  } else {
    paste("WGCNA's pickSoftThreshold recommends a power threshold of ",
          sft()$powerEstimate,
          "which corresponds to a R squared for SFT of ",
          round(sft()$fitIndices[sft()$fitIndices$Power ==
                                   sft()$powerEstimate, "SFT.R.sq"],2))
  }
})

# Display only if the table is displayed
output$sft_table_title <- renderUI({
  req(sft())
  HTML(paste0("<h2> WGCNA's pickSoftThreshold fitIndices table </h2>"))
})


output$sft_table <- renderTable({
  req(sft())
  return(sft()$fitIndices %>% as.data.frame())},
  striped = TRUE,
  bordered = TRUE)


# Creates a correspondance table with modules pre and postmerging
modules_pre_post <- reactive({
  req(modules())
  module_post_df <- stack(modules()$modules)
  
  module_pre_df <- stack(modules()$modules_premerge)
  
  gr_joined <- dplyr::left_join(module_pre_df, module_post_df, by = "values") %>%
    dplyr::group_by(ind.y) %>%
    summarise(modules_pre = unique(ind.x))
  
  colnames(gr_joined) <- c("Modules_post", "Modules_pre")
  aggregate(Modules_pre~Modules_post, gr_joined, paste, collapse=", ")
})


output$merge_premerge <- renderTable({
  modules_pre_post()
}, bordered = TRUE,
striped = TRUE)
  

output$merge <- renderPlot({
  req(modules())
 plot_modules_merge(
    modules_premerge = modules()$modules_premerge, 
    modules_merged = modules()$modules,
    zoom = 1)
  
}, res = 96, height = 600)


output$sizes <- renderPlot({
ggplot(data.frame(modules()$modules %>% stack), 
                aes(x = ind)) + stat_count() +
  ylab("Number of genes") +
  xlab("Module")
}, res = 96)

observe({
  updateSelectizeInput(session,
                       inputId = "select_mod", 
                       choices = names(modules()$modules),
                       # To select the first module
                       selected = names(modules()$modules)[2],
                       server = TRUE)
})


output$Enrichment <- renderPlotly({
  req(modules_enriched())
  # Plot only if there is something to plot
  if(input$select_mod %in% modules_enriched()$result$query) {
    plot_enrichment(modules_enriched(),
                    modules = input$select_mod)
  }
})

# to be able to write to csv one needs to collapse a list
enrich_down <- reactive({
  tmp <- modules_enriched()$result %>% filter(query == input$select_mod)
  tmp$parents <- vapply(tmp$parents, paste, collapse = ", ", character(1L))
  tmp
})

output$download_enr <- downloadHandler(
  filename = function() {
    paste("enrichment_module_", input$select_mod, ".csv", sep = "")
  },
  content = function(file) {
    write.csv(enrich_down(), file)
  }
)