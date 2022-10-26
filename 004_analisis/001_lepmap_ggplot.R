

######################
### ggplot LepMap3 ###
######################


# ñlean up of the memory
rm(list = ls())

##library
library(tidyverse)

## create a variable of main directory
main_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

# creating a folder of boxplot results
dir.create(file.path(main_dir, "/done/accuracy/"))

##
x11()

## create a list of chromosomes
num_LG <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
            "11", "12", "13", "14", "15", "16", "17", "18", "19", "20") #change

## loop to process each chromosome
for (LG in num_LG) {

  ## import data
  setwd(paste(main_dir, "/data/mapping/dixie_on_noble_LG", LG, "/", sep="")) #change
  order <- read.csv("order.txt", sep = "", skip = 2)
  cut_call <- read.csv("cut_call.txt", sep = "")
  call <- cut_call[7:nrow(cut_call),]

  ## empty variables for loop
  j <- c()
  combined <- c()
  accuracy <- c()

  ## cycle of creating a date frame for plot
  for (i in 1:nrow(order)) {
    j <- row.names(order[i, ])
    j <- as.numeric(j)
    combined <- cbind(call[j, 2], order[i, 2])
    accuracy <- rbind(accuracy, combined)
  }

  ## chromosome flip
  if (as.numeric(accuracy[10, 1]) > as.numeric(accuracy[nrow(accuracy), 1])) {
    roll_90 <- t(apply(accuracy, 2, rev))
    roll_180 <- t(roll_90)
    vector_roll <- as.numeric(roll_180[, 2])
    for (i in 1:nrow(accuracy)) {
      roll_180[i, 2] <- vector_roll[1] - vector_roll[i]
    }
    accuracy <- roll_180
  }

  ## variable for plot
  physical <- as.numeric(accuracy[, 1])
  genetic <- as.numeric(accuracy[, 2])

  ## creating plot
  ggplot(data = NULL, aes(x = physical, y = genetic)) + geom_point()
  ggsave(paste(main_dir, "/done/accuracy/dixie_on_noble_LG", LG, ".tiff", sep = "")) #change

}

setwd(main_dir)
