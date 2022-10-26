

######################
### ggplot LepMap3 ###
######################


# ñlean up of the memory
rm(list = ls())

##library
library(tidyverse)

## create a variable of main directory
main_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

## import data
setwd(paste(main_dir, "/done/quantitative/", sep="")) #change
order <- read.csv("quantitative_position.csv", sep = ";",)

## create a list of chromosomes
num_LG <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
            "11", "12", "13", "14", "15", "16", "17", "18", "19", "20") #change

## loop to process each chromosome
for (LG in num_LG) {
  order[, 1] <- str_replace(order[, 1], paste("dixie_LG", LG, "_marker_", sep = ""), "")
}

order[, 3] <- gsub(",", ".", order[, 3])

##
#order[, 1] <- as.numeric(order[, 1]) / 1000000

cor <- c()

## loop to process each chromosome
for (LG in c(1:20)) {
  
  accuracy <- subset(order, order[, 2] == LG)
  Mb <- as.numeric(accuracy[, 1])
  cM <- as.numeric(accuracy[, 3])

  #ggplot(data = NULL, aes(x = Mb, y = cM)) + geom_point() + ggtitle(paste("Chr", LG, sep = " "))

  a <- cor.test(cM, Mb, method = "pearson")

  cor <- cbind(cor, a$estimate)

}

setwd(main_dir)
