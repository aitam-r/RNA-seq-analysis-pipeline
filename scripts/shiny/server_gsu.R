coldata <- reactive({
  req(input$upload_samp_tab)
  
  extension <- tools::file_ext(input$upload_samp_tab$name)
  tmp <- switch(extension,
         csv = vroom::vroom(input$upload_samp_tab$datapath, delim = ","),
         tsv = vroom::vroom(input$upload_samp_tab$datapath, delim = "\t"),
         validate("Invalid file : Need a .csv or .tsv file")
  )
  validate(need("names" %in% colnames(tmp), "No column named names"))
  tmp
})

observe({
  req(input$upload_counts)
  
  extension <- tools::file_ext(input$upload_counts$name)
  my_values$counts <- switch(extension,
                             csv = vroom::vroom(input$upload_counts$datapath, delim = ","),
                             tsv = vroom::vroom(input$upload_counts$datapath, delim = "\t"),
                             validate("Invalid file : Need a .csv or .tsv file")
  )
  
})

observe({
  # make gene ids rownames
  if(!tibble::has_rownames(req(my_values$counts)))
    my_values$counts <- tibble::column_to_rownames(my_values$counts, colnames(my_values$counts)[1])
}) 