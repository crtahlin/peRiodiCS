# prepare viral data for usage - import from txt file to .rda file and put into /data folder 
viral_east_mediteranean <- read.delim("data-raw/viral_east_mediteranean.txt", sep="\t", 
             na.strings = c(999))

devtools::use_data(viral_east_mediteranean, overwrite = TRUE)
