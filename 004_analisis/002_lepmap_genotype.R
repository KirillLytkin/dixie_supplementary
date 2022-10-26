

########################
### genotype LepMap3 ###
########################


# ñlean up of the memory
rm(list = ls())

##library
library(stringr)

## create a variable of main directory
main_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

# creating a folder of dones
dir.create(file.path(main_dir, "/done/genotype/"))

## create a list of chromosomes
num_LG <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
            "11", "12", "13", "14", "15", "16", "17", "18", "19", "20") #change

## input phenotype
setwd(paste(main_dir, "/data/phenotype/", sep=""))
data_phenotype <- read.csv("phenotype.txt", sep = "")

##
for (phe in 1:nrow(data_phenotype)){
  
  #phe <- 1
  
  ## empty variable
  All_LG_QTL <- c()
  
  ##
  phenotype_name <- data_phenotype[phe, 1]
  dir.create(file.path(main_dir, "/done/genotype/", phenotype_name))
  phenotype <- data_phenotype[phe, 2:ncol(data_phenotype)]
  
  ## loop to process each chromosome
  for (LG in num_LG) {
    
    #LG <- "02"
    
    ## import data
    setwd(paste(main_dir, "/data/mapping/dixie_on_noble_LG", LG, "/", sep="")) #change
    data_genotype_code <- read.csv("genotype.txt", sep = "", skip = 5)
    
    data_genotype <- data_genotype_code[7:ncol(data_genotype_code)]
    #data_cM <- data_genotype_code[4]
    
    setwd(paste(main_dir, "/data/mapping/", sep=""))
    data_pedigree <- read.csv("pedigree.txt", sep = "", header = FALSE)
    
    ## for physical position
    setwd(paste(main_dir, "/data/mapping/dixie_on_noble_LG", LG, "/", sep="")) #change
    order <- read.csv("order.txt", sep = "", skip = 2)
    cut_call <- read.csv("cut_call.txt", sep = "")
    call <- cut_call[7:nrow(cut_call),]
    
    ## cut genotype
    data_genotype <- data_genotype[1:nrow(order),]
    
    ## empty variables for loop
    genotype <- c()
    
    ## loop to decode
    for (i in 1:nrow(data_genotype)) {
      a <- str_replace(data_genotype[i, ], "11", "a")
      b <- str_replace(a, "22", "a")
      h1 <- str_replace(b, "12", "h")
      h2 <- str_replace(h1, "21", "h")
      genotype <- rbind(genotype, h2)
    }
    
    ## data generation for rQTL
    first_row <- cbind("id", "V2", "V3", data_pedigree[2, 5:ncol(data_pedigree)])
    quantitative <- cbind(phenotype_name, "V2", "V3", phenotype)
    colnames(quantitative) <- first_row
    marker <- paste("dixie_LG", LG, "_marker_", sep = "") #change
    
    ## empty variables for loop
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
      
      roll_90 <- t(apply(genotype, 2, rev))
      roll_180 <- t(roll_90)
      genotype <- roll_180
    }
    
    nrow(accuracy)
    nrow(genotype)
    
    ## loop to create data for rQTL
    for (i in 1:nrow(genotype)) {
      rest_row <- c(paste(marker, accuracy[i, 1], sep = ""), as.numeric(LG), accuracy[i, 2])
      rest_row <- t(rest_row)
      for (j in 1:ncol(genotype)) {
        rest_row <- cbind(rest_row, genotype[i, j])
      }
      colnames(rest_row) <- first_row
      quantitative <- rbind(quantitative, rest_row)
    }
    
    ## transform and prepare data
    quantitative <- t(quantitative)
    
    ## renaming col
    colnames(quantitative) <- quantitative[1, ] 
    quantitative <- quantitative[2:nrow(quantitative), ]
    
    ## first col for All_LG_QTL
    first_col <- quantitative[, 1]
    
    ## empty case
    quantitative[1:2, 1] <- ""
    
    ## export of received data to done (if needed)
    #setwd(paste(main_dir, "/done/genotype/", sep=""))
    #write.csv(quantitative,
    #paste("dixie_on_noble_LG", LG, "_lepmap.csv", sep=""),
    #row.names = FALSE) #change
    
    ## combining chromosomes into one file
    All_LG_QTL <- cbind(All_LG_QTL, quantitative[, 2:ncol(quantitative)])
    
  }
  
  ## adding the first col
  All_LG_QTL <- cbind(first_col, All_LG_QTL)
  colnames(All_LG_QTL)[1] <- colnames(quantitative)[1]
  
  ## empty case
  All_LG_QTL[1:2, 1] <- ""
  
  ## id col
  id <- c("", "", 1:nrow(All_LG_QTL))
  All_LG_QTL <- cbind(id, All_LG_QTL)
  
  ## save dones
  setwd(paste(main_dir, "/done/genotype/", phenotype_name, "/",sep=""))
  write.csv(All_LG_QTL, "quantitative.csv", row.names = FALSE)
  
}

## save dones
visual <- t(All_LG_QTL)

## replace "." to "," 
replace_visual <- gsub("\\.", ",", visual)

setwd(paste(main_dir, "/done/genotype/", sep=""))
write.csv2(replace_visual, "visual.csv")

setwd(main_dir)
