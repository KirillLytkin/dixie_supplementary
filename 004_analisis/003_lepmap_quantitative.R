

###########
### QTL ###
###########


## erase all graphs
graphics.off()

## ñlean up of the memory
rm(list=ls())

## library
library(qtl)

## create a variable of main directory
main_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

## creating a folder of quantitative
dir.create(file.path(main_dir, "/done/quantitative/"))

## input phenotype for loop
setwd(paste(main_dir, "/data/phenotype/", sep=""))
data_phenotype <- read.csv("phenotype.txt", sep = "")

##
x11()

##
for (phe in 1:nrow(data_phenotype)){
  
  #phe <- 1
  
  ## creating a folder of done
  phenotype_name <- data_phenotype[phe, 1]
  dir.create(file.path(main_dir, "/done/quantitative/", phenotype_name))

  ## import data
  setwd(paste(main_dir, "/done/genotype/", phenotype_name, "/", sep=""))

  ## load the .csv file with the genotype and phenotype data
  phenoGeno_0 <- read.cross(format="csv",
                            file="quantitative.csv",  
                            na.strings=c("NA"),
                            genotypes=c("a", "h"))
  ## warning message: some markers at the same position on chr
  ## 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20; use jittermap().
  
  ## use jittermap
  phenoGeno <- jittermap(phenoGeno_0, amount=1e-6)
  
  ## recalculation of genetic positions
  newmap <- est.map(phenoGeno)
  logliks <- sapply(newmap, attr, "loglik")
  
  ## building and saving a plot
  setwd(paste(main_dir, "/done/quantitative/", phenotype_name, "/", sep=""))
  plotMap(phenoGeno, newmap)
  dev.copy(tiff, filename="comparison.tiff") 
  dev.off()
  
  ## replacement of genetic positions
  phenoGeno <- replace.map(phenoGeno, newmap)

  ## summary of the data file and plotting
  #summary(phenoGeno)
  capture.output(summary(phenoGeno), file = "summary.txt")
  
  plot(phenoGeno)
  dev.copy(tiff, filename="summary.tiff") 
  dev.off()
  
  ## calculating genotype probabilities
  phenoGeno <- calc.genoprob(phenoGeno,
                             step=0,
                             off.end=0.0,
                             error.prob=1.0e-4,
                             stepwidth="fixed",
                             map.function="kosambi")
  
  phenoGeno_cross <- sim.geno(phenoGeno,
                              n.draws=32,
                              step=0,
                              off.end=0.0,
                              error.prob=1.0e-4,
                              stepwidth="fixed",
                              map.function="kosambi")

  ## running QTL scan using CIM method with permutations
  scan.cim <- cim(phenoGeno_cross,
                  pheno.col=2,
                  map.function="kosambi")
  ## warning messages: X'X matrix is singular.
  
  ## type the below command in the console to run the qtl scan using CIM method
  scan.cim.perm <- cim(phenoGeno_cross,
                       pheno.col=2,
                       map.function="kosambi",
                       n.perm=1000)
  ## warning messages: X'X matrix is singular.
  
  ## threshold
  value <- 0.05 #change
  
  ## run permutation test by typing the below command
  sum.scan.sim.perm <- summary(scan.cim.perm)
  sum.scan.sim <- summary(scan.cim, threshold=value)
  
  marker <- subset(sum.scan.sim, sum.scan.sim$lod == max(sum.scan.sim$lod))
  #marker.perm <- subset(sum.scan.sim, sum.scan.sim$lod > sum.scan.sim.perm[1])
  lod.perm <- round(sum.scan.sim.perm[1], 2)
  
  ## plot the QTL scan
  plot(scan.cim)
  abline(h=lod.perm, col="blue")
  text(2800, round(marker[3], 1), paste("lod.perm = ", lod.perm, sep = ""))
  dev.copy(tiff, filename="scan.tiff") 
  dev.off()
  
  ##
  scan.cim.chr <- c()
  
  for (i in 1:nrow(scan.cim)) {
    if (scan.cim[i, 1] == as.numeric(marker[1])) {
      scan.cim.chr <- rbind(scan.cim.chr, scan.cim[i,])
    }
  }
  
  ## plot the QTL scan chr
  plot(scan.cim.chr)
  abline(h=lod.perm, col="blue")
  title(paste("Ñhromosome", marker[1], "containing a significant QTL"))
  text(140, round(marker[3], 1), paste("lod.perm = ", lod.perm, sep = ""))
  dev.copy(tiff, filename="scan chromosome.tiff") 
  dev.off()
  
  ## QTL effect plot
  plotPXG(phenoGeno_cross,
          pheno.col=2,
          marker=row.names(marker))
  dev.copy(tiff, filename="effect.tiff") 
  dev.off()
  
  ## make QTL model with the significant marker
  qtl <- makeqtl(phenoGeno_cross,
                 chr=marker[1],
                 pos=marker[2],
                 what=c("prob")) 

  fitqtl <- fitqtl(phenoGeno_cross,
                   pheno.col=2,
                   formula=y~Q1,
                   qtl=qtl,
                   method="hk",
                   get.ests=T)
  ## warning message: Dropping 16 individuals with missing phenotypes.
  
  ##
  #sum.fitqtl <- summary(fitqtl)
  capture.output(summary(fitqtl), file = "fitqtl.txt")
  
  ## identify QTL intervals using LOD drop
  scan.lodint <- lodint(scan.cim, chr=as.numeric(marker[1]), drop=1.8)
  capture.output(scan.lodint, file = "lodint.txt")
  
}


## save position for plot
setwd(paste(main_dir, "/done/quantitative/", sep=""))
write.csv2(scan.cim, "quantitative_position.csv")

##
setwd(main_dir)
