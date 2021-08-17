server <- function(input, output, session) {
  
  # Object to store reactives that can change
  my_values <- reactiveValues(coldata = NULL,
                              se = NULL,
                              variable_chosen = NULL,
                              counts = NULL,
                              res = NULL,
                              counts_norm = NULL,
                              counts_filt = NULL)

  source(file = "server_gsu.R", local = TRUE)
  # add an object to check user input
  iv <- InputValidator$new()
  source(file = "server_deseq2.R", local = TRUE)
  source(file = "server_wgcna.R", local = TRUE)  
}
