
rm(list = ls())

library(ggplot2)
library(dplyr)

main_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

setwd(paste(main_dir, "/data/", sep=""))
df <- read.csv("dixie_on_noble.recode.hmp.txt", sep = "")
df_after <- read.csv("after_filtration.txt", sep = "") 

LG_num <- c("CP092926.1", "CP092927.1", "CP092928.1", "CP092929.1", "CP092930.1",
            "CP092931.1", "CP092932.1", "CP092933.1", "CP092934.1", "CP092935.1",
            "CP092936.1", "CP092937.1", "CP092938.1", "CP092939.1", "CP092940.1",
            "CP092941.1", "CP092942.1", "CP092943.1", "CP092944.1", "CP092945.1")

setwd(paste(main_dir, "/done/", sep=""))
tiff("alllg_follow.tiff", units="in", width=10, height=8, res=300)
par(mfrow = c(4, 5))

for (i in 1:20) {
  
  chrom_pos <- c()
  chrom_pos <- cbind(df$chrom, df$pos)
  
  LG <- LG_num[i]
  
  page <- subset(chrom_pos, chrom_pos[, 1] == LG)
  page_after <-  subset(df_after, df_after[, 1] == LG)
  
  level <- rep(1.2, nrow(page))
  level_after <- rep(1.8, nrow(page_after))
  x_zero <- c(0, 20)
  y_zero <- rep(1.5, 2)
  
  page <- cbind(page, level)
  page_after <- cbind(page_after, level_after)
  
  page[, 2] <- as.numeric(page[, 2]) / 1000000
  page_after[, 2] <- as.numeric(page_after[, 2]) / 1000000
  
  plot(x_zero, y_zero, col = "white", frame = FALSE, yaxt='n',
       xlab = "Mb", ylab = "before | after", main = paste("LG", i))
  lines(page[, 2], page[, 3], pch = "|", type = "o", col = "blue")  
  lines(page_after[, 2], page_after[, 3], pch = "|", type = "o", col = "blue")
  
}

dev.off()
setwd(main_dir)
