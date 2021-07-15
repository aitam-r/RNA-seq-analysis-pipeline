library(tximeta)

setwd("~/Documents/TAAAAF/Stage/liver/scripts/shiny/")
head_dir <- "../.."
dir_quants <- file.path(head_dir, "data/quants")
coldata <- read.table(file.path(head_dir, "coldata.csv"), header = T, sep = ",")
coldata$files <- file.path(dir_quants, coldata$names, "quant.sf")
stopifnot(file.exists(coldata$files))

##  Annotation et importation
se <- tximeta(coldata)

## Only gene expression is of interest here
gse <- summarizeToGene(se)

## Arrangement des facteurs
gse$condition %<>% factor()




save.image("./txidata.RData")