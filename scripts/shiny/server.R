server <- function(input, output, session) {
  
  # add an object to check user keyboard input
  iv <- InputValidator$new()
  source(file = "server_deseq2.R", local = TRUE)
  source(file = "server_wgcna.R", local = TRUE)  
}
