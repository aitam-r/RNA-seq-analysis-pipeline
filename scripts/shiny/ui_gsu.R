tab_gsu <- tabItem("gsu",
                   fluidRow(
                     column(width = 6,
                       checkboxInput(inputId = "snakemake",
                                     label = "I have used the built-in RNA quantification",
                                     value = TRUE),
                       br(),
                       p("Please provide the sample data table (.csv or .tsv). It should have a column containing identifiers (called names) for each sample, and then variable(s) in further columns"),
                       p("Furthermore, the column names should not have spaces or -"),
                       br(),
                       fileInput("upload_samp_tab", NULL, accept = c(".csv", ".tsv"))
                     ),
                     column(width = 6,
                            conditionalPanel(condition = "input.snakemake == false",
                                             br(),
                                             p("Please provide the count tables, column names being the sample names. They should imperatively match the names in the sample data table. They should also be as close as possible to raw counts."),
                                             br(),
                                             fileInput("upload_counts", NULL, accept = c(".csv", ".tsv"))
                            ),
                            conditionalPanel(condition = "input.snakemake == true",
                                             actionButton("import",
                                                          label = "Import Salmon count data")
                            )
                     )
                   )
)