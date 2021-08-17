observeEvent(input$upload_samp_tab, {
  req(input$upload_samp_tab)
  
  extension <- tools::file_ext(input$upload_samp_tab$name)
  my_values$coldata <- switch(extension,
         csv = vroom::vroom(input$upload_samp_tab$datapath, delim = ","),
         tsv = vroom::vroom(input$upload_samp_tab$datapath, delim = "\t"),
         validate("Invalid file : Need a .csv or .tsv file")
  )
  validate(need("names" %in% colnames(my_values$coldata), "No column named names"))
  my_values$coldata %<>% as.data.frame()
  for (i in 2:ncol(my_values$coldata)) {
    my_values$coldata[, i] %<>% as.factor()
  }
  # my_values$coldata %<>% tibble::column_to_rownames("names")
})

observeEvent(input$upload_counts, {
  
  extension <- tools::file_ext(input$upload_counts$name)
  my_values$counts <- switch(extension,
                             csv = vroom::vroom(input$upload_counts$datapath, delim = ","),
                             tsv = vroom::vroom(input$upload_counts$datapath, delim = "\t"),
                             validate("Invalid file : Need a .csv or .tsv file")
  )
  # make gene ids rownames
  if(!tibble::has_rownames(req(my_values$counts)))
    my_values$counts <- tibble::column_to_rownames(my_values$counts, var = colnames(my_values$counts)[1])
})

observeEvent(input$import, {
  req(my_values$coldata)
  # head_dir <- "../.."
  # dir_quants <- file.path(head_dir, "data/quants")
  my_values$coldata$files <- file.path(getwd(),
                                       "../../data/quants",
                                       my_values$coldata$names,
                                       "quant.sf")
  validate(need(file.exists(my_values$coldata$files),
               message = "Some files specified by sample data names do not exist"))
  if(!txi_met_chosen) {
    withProgress(message = "Importing Salmon counts, GENCODE", {
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
      keys <- keys(txdb, keytype = "TXNAME")
      tx2gene <- AnnotationDbi::select(txdb, keys, "GENEID", "TXNAME")
      txi <- tximport(my_values$coldata$files,
                      type = "salmon",
                      tx2gene = tx2gene)
      validate(need(all(my_values$coldata$names == colnames(txi$counts))),
               "Names in sample data table do not correspond to names of salmon quants")
    })
  }
  else {
    withProgress(message = "Importing Salmon counts, tximeta/Ensembl", {
    my_values$se <- tximeta(my_values$coldata)
    })
  }
})

gse <- eventReactive(my_values$se, {
  withProgress(message = "Summarising to gene", {
  gse <- summarizeToGene(my_values$se)
  gse$condition %<>% factor()
  })
  return(gse)
})