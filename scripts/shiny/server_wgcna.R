# Calculs ------

my_values$counts_norm <- reactive({
  #using bias corrected counts (without an offset)
  raw_counts <- assay(gse_scaled, "counts") %>% round()
  # Remove non-expressed genes
  raw_counts <- raw_counts[rowSums(raw_counts) > 0,]
  
  # Calculate upper quartile value (not taking into account 0s by sample)
  quantile_expressed <- apply(raw_counts, 
                              2, 
                              function(x){quantile(x[x>0], 0.75)})
  # Divide each column by its upper quartile value
  t(raw_counts)/quantile_expressed
})

observeEvent(input$rm_sample, {
  my_values$counts_norm() <- 
    my_values$counts_norm()[!(rownames(my_values$counts_norm()) %in%
                              input$rm_sample),]
})


net <- eventReactive(input$build, {
  withProgress(message = "Building network", {
    build_net(my_values$counts_filt(),
              cor_func = "spearman",
              power_value = input$sft_thres,
              network_type = isolate(input$type_net),
              n_threads = 1)
  })
})
modules <- eventReactive(net(), {
  withProgress(message = "Detecting Modules", {
    detect_modules(my_values$counts_filt(),
                   net()$network,
                   detailled_result = T)
  })
})

modules_enriched <- eventReactive({
  modules()
  input$sources
  input$enrich}, { 
    withProgress(message = "Performing Enrichment", {
      bio_enrich(modules()$modules,
                 organism = "hsapiens",
                 sources = input$sources)
    })
  })

observeEvent(my_values$counts_norm(), {
  # Update the maximum of gene one can keep
  # default value of 8 000 genes not to crush computer resources
  m <- round(8000/ncol(req(my_values$counts_norm())), 1)
  updateSliderInput(session,
                    "percent_g",
                    max = m)
})

my_values$counts_filt <- eventReactive(my_values$counts_norm(), {
  filter_low_var(my_values$counts_norm(),
                 pct = input$percent_g,
                 type = "median")
})

output$outliers <- renderPlot({
  req(input$explore_w)
  sampleTree <- hclust(dist(my_values$counts_norm()), method = "average")
  plot(sampleTree, main = "Sample clustering to detect outliers",
       sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
})

output$warning <- renderText({
  req(my_values$counts_norm())
  if(nrow(my_values$counts_norm()) < 20) {
    paste("\n", "Warning : The number of samples is too low for WGCNA analyis.", "\n")
  }
})

sft <- eventReactive({input$update_sft}, {
  withProgress(message = "Calculating", {
    pwr_vec <- c(1:9, seq(10, 30, by = 2))
    WGCNA::pickSoftThreshold(req(my_values$counts_filt()),
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

output$sft_table_title <- renderUI({
  req(sft())
  HTML(paste0("<h2> WGCNA's pickSoftThreshold fitIndices table </h2>"))
})


output$sft_table <- renderTable({
  req(sft())
  return(sft()$fitIndices %>% as.data.frame())},
  striped = TRUE,
  bordered = TRUE)


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
