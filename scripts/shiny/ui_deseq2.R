tab_deseq2_su <- tabItem("deseq2_su",
                         #Choose the experimental design
                         selectInput(inputId = "variables",
                                     label = "Choose the variables of the experimental design :",
                                     choices = NULL,
                                     multiple = TRUE),
                        
                         # Choose the variable across which differential expression is calculated 
                         selectInput(inputId = "deseq_var",
                                     label = "Choose the variable for differential Expression",
                                     choices = NULL),
                                     
                         
                         # Choose base level for condition
                         selectInput(inputId = "base_cond",
                                     label = "Choose the base condition :",
                                     choices = NULL),
                         
                         # Choose base level for condition
                         selectInput(inputId = "compare_cond",
                                     label = "Choose the condition to compare with:",
                                     choices = NULL,
                                     selected = NULL),
                         
                         # rld or vst?
                         selectInput(inputId = "r_v",
                                     label = "Choose the variance stabilizing transformation wanted for the PCA (vst is much faster) :",
                                     choices = c("rlog", "vst"),
                                     selected = "vst"),
                         
                         actionButton(inputId = "execute_d",
                                      label = "Run DESeq2")
)
                        

tab_plots <- tabItem("plots",
                     tabBox(width = NULL,
                       tabPanel("Sample-to-sample distances", 
                                plotOutput(outputId = "dist")),
                       
                       tabPanel("PCA", 
                                plotOutput(outputId = "pca")),
                       
                       tabPanel("MAplot", 
                                plotOutput(outputId = "ma")),
                       
                       tabPanel("Volcano Plot", 
                                plotOutput(outputId = "volcano"))
                       )
)
                        

tab_deg <- tabItem("deg",
                   sidebarLayout(
                     
                     sidebarPanel(
                       numericInput(inputId = "pval_cutoff",
                                    label = "Enter the maximum p-value :",
                                    value = 0.05,
                                    min = 0,
                                    max = 1,
                                    step = .05),
                       
                       numericInput(inputId = "lfc_cutoff",
                                    label = "Enter the minimum (absolute) logFoldChange :",
                                    value = 1,
                                    min = 0)
                     ),
                     
                     mainPanel(
                       
                       textOutput(outputId = "nb_genes"),
                       br(),
                       htmlOutput("genes_na"),
                       br(),
                       DT::dataTableOutput(outputId = "genes"),
                       
                       downloadButton(outputId = "download",
                                      label = "Download table")
                     )
                   )
)
                        

tab_sgc <- tabItem("sgc",
                   sidebarLayout(
                     sidebarPanel(
                       selectizeInput(inputId = "sel_gene",
                                      label = "Select which gene to plot :",
                                      choices = NULL),
                       
                       checkboxGroupInput("condition_plot",
                                          label = "Choose which condition to plot :",
                                          choices = NULL,
                                          selected = NULL),
                       
                       selectInput("plot",
                                   label = "Select the type of plot",
                                   choices = c("base", "barplot", "boxplot"),
                                   selected = "base"),
                       
                       sliderInput("alpha",
                                   label = "Choose the barplot transparency",
                                   min = 0,
                                   max = 1,
                                   value = .2)
                     ),
                     mainPanel(
                       plotOutput(outputId = "plot_gene"),
                       
                       downloadButton(outputId = "down_plot",
                                      label = "Download Plot")
                     )
                   )
)
