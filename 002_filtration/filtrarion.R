

#########################
### parent filtration ###
#########################


## create a variable of main directory
main_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

## import data
setwd(paste(main_dir, "/data/", sep=""))
df <- read.csv("dixie_on_noble.recode.hmp.txt", sep = "")

## address the area with data of parents
df_parent <- df[14:15]

## count the number of each heterozygous
AC <- rowSums(df_parent == "AC")
AG <- rowSums(df_parent == "AG")
AT <- rowSums(df_parent == "AT")
CA <- rowSums(df_parent == "CA")
CG <- rowSums(df_parent == "CG")
CT <- rowSums(df_parent == "CT")
TA <- rowSums(df_parent == "TA")
TG <- rowSums(df_parent == "TG")
TC <- rowSums(df_parent == "TC")
GA <- rowSums(df_parent == "GA")
GC <- rowSums(df_parent == "GC")
GT <- rowSums(df_parent == "GT")
df_parent_hete <- cbind(AC, AG, AT,
                        CA, CG, CT,
                        TA, TG, TC,
                        GA, GC, GT)
df_parent_hete_sum <- rowSums(df_parent_hete)

df_parent_filtered <- subset(df, df_parent_hete_sum >= 1)


#############
### chisq ###
#############


## address the area with data of individuals
df_ind <- df_parent_filtered[20:105]

## count the number of each heterozygous
AA <- rowSums(df_ind == "AA")
GG <- rowSums(df_ind == "GG")
TT <- rowSums(df_ind == "TT")
CC <- rowSums(df_ind == "CC")

## count the number of each heterozygous
AC <- rowSums(df_ind == "AC")
AG <- rowSums(df_ind == "AG")
AT <- rowSums(df_ind == "AT")
CA <- rowSums(df_ind == "CA")
CG <- rowSums(df_ind == "CG")
CT <- rowSums(df_ind == "CT")
TA <- rowSums(df_ind == "TA")
TG <- rowSums(df_ind == "TG")
TC <- rowSums(df_ind == "TC")
GA <- rowSums(df_ind == "GA")
GC <- rowSums(df_ind == "GC")
GT <- rowSums(df_ind == "GT")

## combine the values of hete- and homozygous
hete_homo <- cbind(AC, AG, AT, CA,
                   CG, CT, TA, TG,
                   TC, GA, GC, GT,
                   AA, GG, TT, CC)

## create a variable to store result
df_filtered <- c()

## calculate chisq in a loop
for (i in 1:nrow(hete_homo)){
  
  ## calculation if 1 hete- and 2 homozygous, i.e. 2 to 1 to 1 segregation
  if (sum(hete_homo[i, 1:12] == 0) == 11 & sum(hete_homo[i, 13:16] == 0) == 2){
    pos_val <- subset(hete_homo[i,], hete_homo[i,] > 0)
    chisq_list <- chisq.test(pos_val, p=c(2/4, 1/4, 1/4))
  
  if (chisq_list$p.value >= 0.05){
    df_filtered <- rbind(df_filtered, df_parent_filtered[i,])
    }
  }
}


##############
### lepmap ###
##############


parent <- df_filtered[14:15]

for (i in 1:nrow(parent)) {
  if (parent[i, 1] == "TT"){
    parent[i, 1] <- parent[i, 2]
  }
  if (parent[i, 1] == "GG"){
    parent[i, 1] <- parent[i, 2]
  }
  if (parent[i, 1] == "CC"){
    parent[i, 1] <- parent[i, 2]
  }
  if (parent[i, 1] == "AA"){
    parent[i, 1] <- parent[i, 2]
  }
  if (parent[i, 1] == "NN"){
    parent[i, 1] <- parent[i, 2]
  }
}

for (i in 1:nrow(parent)) {
  if (parent[i, 2] == "TT"){
    parent[i, 2] <- parent[i, 1]
  }
  if (parent[i, 2] == "GG"){
    parent[i, 2] <- parent[i, 1]
  }
  if (parent[i, 2] == "CC"){
    parent[i, 2] <- parent[i, 1]
  }
  if (parent[i, 2] == "AA"){
    parent[i, 2] <- parent[i, 1]
  }
  if (parent[i, 2] == "NN"){
    parent[i, 2] <- parent[i, 1]
  }
}

df_lepmap <- df_filtered
df_lepmap[14:15] <- parent

library(stringr)
colnames(df_lepmap) <- str_replace(colnames(df_lepmap),
                                   ".trimmed.fq.gz", "")


##########
### LG ###
##########


## export of received data to result
setwd(paste(main_dir, "/done/", sep=""))
write.csv2(subset(df_lepmap[,c(1:11,14:15,20:105)],
                  df_lepmap$chrom == "CP092926.1"),
           "dixie_on_noble_LG01_lepmap.csv", row.names = FALSE)
write.csv2(subset(df_lepmap[,c(1:11,14:15,20:105)],
                  df_lepmap$chrom == "CP092927.1"),
           "dixie_on_noble_LG02_lepmap.csv", row.names = FALSE)
write.csv2(subset(df_lepmap[,c(1:11,14:15,20:105)],
                  df_lepmap$chrom == "CP092928.1"),
           "dixie_on_noble_LG03_lepmap.csv", row.names = FALSE)
write.csv2(subset(df_lepmap[,c(1:11,14:15,20:105)],
                  df_lepmap$chrom == "CP092929.1"),
           "dixie_on_noble_LG04_lepmap.csv", row.names = FALSE)
write.csv2(subset(df_lepmap[,c(1:11,14:15,20:105)],
                  df_lepmap$chrom == "CP092930.1"),
           "dixie_on_noble_LG05_lepmap.csv", row.names = FALSE)
write.csv2(subset(df_lepmap[,c(1:11,14:15,20:105)],
                  df_lepmap$chrom == "CP092931.1"),
           "dixie_on_noble_LG06_lepmap.csv", row.names = FALSE)
write.csv2(subset(df_lepmap[,c(1:11,14:15,20:105)],
                  df_lepmap$chrom == "CP092932.1"),
           "dixie_on_noble_LG07_lepmap.csv", row.names = FALSE)
write.csv2(subset(df_lepmap[,c(1:11,14:15,20:105)],
                  df_lepmap$chrom == "CP092933.1"),
           "dixie_on_noble_LG08_lepmap.csv", row.names = FALSE)
write.csv2(subset(df_lepmap[,c(1:11,14:15,20:105)],
                  df_lepmap$chrom == "CP092934.1"),
           "dixie_on_noble_LG09_lepmap.csv", row.names = FALSE)
write.csv2(subset(df_lepmap[,c(1:11,14:15,20:105)],
                  df_lepmap$chrom == "CP092935.1"),
           "dixie_on_noble_LG10_lepmap.csv", row.names = FALSE)
write.csv2(subset(df_lepmap[,c(1:11,14:15,20:105)],
                  df_lepmap$chrom == "CP092936.1"),
           "dixie_on_noble_LG11_lepmap.csv", row.names = FALSE)
write.csv2(subset(df_lepmap[,c(1:11,14:15,20:105)],
                  df_lepmap$chrom == "CP092937.1"),
           "dixie_on_noble_LG12_lepmap.csv", row.names = FALSE)
write.csv2(subset(df_lepmap[,c(1:11,14:15,20:105)],
                  df_lepmap$chrom == "CP092938.1"),
           "dixie_on_noble_LG13_lepmap.csv", row.names = FALSE)
write.csv2(subset(df_lepmap[,c(1:11,14:15,20:105)],
                  df_lepmap$chrom == "CP092939.1"),
           "dixie_on_noble_LG14_lepmap.csv", row.names = FALSE)
write.csv2(subset(df_lepmap[,c(1:11,14:15,20:105)],
                  df_lepmap$chrom == "CP092940.1"),
           "dixie_on_noble_LG15_lepmap.csv", row.names = FALSE)
write.csv2(subset(df_lepmap[,c(1:11,14:15,20:105)],
                  df_lepmap$chrom == "CP092941.1"),
           "dixie_on_noble_LG16_lepmap.csv", row.names = FALSE)
write.csv2(subset(df_lepmap[,c(1:11,14:15,20:105)],
                  df_lepmap$chrom == "CP092942.1"),
           "dixie_on_noble_LG17_lepmap.csv", row.names = FALSE)
write.csv2(subset(df_lepmap[,c(1:11,14:15,20:105)],
                  df_lepmap$chrom == "CP092943.1"),
           "dixie_on_noble_LG18_lepmap.csv", row.names = FALSE)
write.csv2(subset(df_lepmap[,c(1:11,14:15,20:105)],
                  df_lepmap$chrom == "CP092944.1"),
           "dixie_on_noble_LG19_lepmap.csv", row.names = FALSE)
write.csv2(subset(df_lepmap[,c(1:11,14:15,20:105)],
                  df_lepmap$chrom == "CP092945.1"),
           "dixie_on_noble_LG20_lepmap.csv", row.names = FALSE)
