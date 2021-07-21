# UI -----------------------------------------------------------------------------

ui <- fluidPage(theme = shinytheme("cosmo"),
                navbarPage(
                  title = "My little app",
                  tabPanel("Setup",
                           #Choose the experimental design
                           selectInput(inputId = "variables",
                                       label = "Choose the variables of the experimental design :",
                                       choices = colnames(coldata),
                                       multiple = T),
                           
                           # Choose base level for condition
                           selectInput(inputId = "base_cond",
                                       label = "Choose the base condition :",
                                       choices = levels(gse$condition)),
                           
                           # Choose base level for condition
                           selectInput(inputId = "compare_cond",
                                       label = "Choose the condition to compare with:",
                                       choices = levels(gse$condition)),
                           
                           actionButton(inputId = "execute",
                                        label = "Run DESeq2")
                  ),
                  
                  tabPanel("Sample-to-sample distances", plotOutput(outputId = "dist")),
                  
                  tabPanel("PCA", plotOutput(outputId = "pca")),
                  
                  tabPanel("MAplot", plotOutput(outputId = "ma")),
                  
                  tabPanel("Volcano Plot", plotOutput(outputId = "volcano")),
                  
                  tabPanel("Most DEGs",
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
                               DT::dataTableOutput(outputId = "genes"),
                               
                               downloadButton(outputId = "download",
                                              label = "Download table")
                             )
                           )
                  ),
                  
                  tabPanel("Plot Gene Count", 
                           sidebarLayout(
                             sidebarPanel(
                               selectizeInput(inputId = "sel_gene",
                                              label = "Select which gene to plot :",
                                              choices = NULL)
                             ),
                             mainPanel(
                               plotOutput(outputId = "plot_gene")
                             )
                           )
                  )
                )
)
