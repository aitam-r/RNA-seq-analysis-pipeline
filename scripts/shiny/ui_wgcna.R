tab_wgcna_su <- tabItem(tabName = "wgcna_su",
                        sidebarLayout(
                          
                          sidebarPanel(
                            br(),
                            # Warning on the number of samples
                            textOutput(outputId = "warning"),
                            br(),
                            actionButton(inputId = "explore_w",
                                         label = "Explore Samples"),
                            br(),br(),
                            
                            selectizeInput(inputId = "rm_sample",
                                           label = "Choose samples to exclude : ",
                                           choices = NULL,
                                           multiple = T,
                                           options = NULL), #list(maxItems = length(coldata$names) - 3)),
                            
                            # Selects the most variable (across all samples) genes. 
                            h6(paste0("The computational load is proportional to square the number of genes. ",
                               "Depending on your computational resources, you might want to limit this number")),
                            sliderInput(inputId = "number_g",
                                        label = "Number of genes used",
                                        value = 8000,
                                        min = 0,
                                        max = 10000,
                                        step = 50),
                            
                            selectInput(inputId = "type_net",
                                        label = "Type of Network",
                                        choices = c("unsigned", "signed", "signed hybrid"),
                                        selected = "signed hybrid"),
                            
                            actionButton(inputId = "update_sft",
                                         label = "Update Power Picked"),
                            
                            br(), br(),br(),
                            
                            numericInput(inputId = "sft_thres", 
                                         label = "Pick soft threshold",
                                         min = 1,
                                         max = 30,
                                         value = 6),
                            
                            textOutput("threshold"),
                            br(),br(),
                        
                            
                            actionButton(inputId = "build",
                                         label = "Build Network")
                            
                          ),
                          mainPanel(
                            plotOutput("outliers"),
                            htmlOutput("sft_table_title"),
                            tableOutput(outputId = "sft_table")
                          )
                        )
)

tab_mod_merge <- tabItem("mod_merge",
                         box(title = "Merging of module",
                             width = 6,
                             plotOutput("merge")
                         ),
                         box(title = "Corresponding modules",
                             width = 6,
                             tableOutput("merge_premerge")
                         )
)


tab_mod_count <- tabItem("mod_count",
                         plotOutput("sizes")
                         
                         # Abandonned : the modules wouldn't correspond to GWENA's input.
                         # ),
                         # tabPanel("TOMplot",
                         #          plotOutput("tomplot")
)


tab_enr <- tabItem("enr",
                   sidebarLayout(
                     sidebarPanel(
                       selectizeInput(inputId = "sources",
                                           label = "Select sources for enrichment",
                                           choices = enr_sources,
                                           selected = c("GO", "KEGG", "REAC"),
                                           multiple = TRUE),
                       
                       actionButton(inputId = "enrich",
                                         label = "Enrich Modules"),
                       br(),br(),
                       
                       selectizeInput(inputId = "select_mod",
                                      label = "Select Module :",
                                      selected = NULL,
                                      choices = NULL)
                     ),
                     mainPanel(
                       plotlyOutput("Enrichment"),
                       
                       downloadButton(outputId = "download_enr",
                                      label = "Download enrichment table")
                     )
                   )
)