library(Rsamtools)
library(tidyverse)
library(data.table)
library(readxl)
library(foreach)
library(doParallel)
# library(foreach)
registerDoParallel(48)

files = list.files(path = paste(getwd(),"mono.poly.UCI.v1.files", sep = "/"), pattern = "*.sam")
files2 = strsplit(files, "_R1_001.sam") %>% sapply(function(x) x[1])
files3 = strsplit(files, "lib3-|_R1_001.sam") %>% sapply(function(x) x[2])

bc.count.df <- data.frame(barcodes = character(0))
a = foreach (i = 1:length(files)) %dopar% {
  Rsamtools::asBam(file = paste(getwd(),"mono.poly.UCI.v1.files",files[i], sep = "/"), destination = paste(getwd(),"mono.poly.UCI.v1.files",files3[i], sep = "/"), overwrite = TRUE)
  bam.df <- as.data.frame(Rsamtools::scanBam(paste0(getwd(),"/","mono.poly.UCI.v1.files","/",files3[i],".bam")))
  bam.df <- bam.df %>% mutate(barcodes = paste(seq,rname,sep = "-"))
  count.dt <- data.table(bam.df[bam.df$flag==0,])[, .N, keyby = bam.df[bam.df$flag==0,]$barcodes]
  colnames(count.dt) <- c("barcodes",files3[i])
  # bc.count.df <- full_join(bc.count.df, count.dt, by = "barcodes")
  return(count.dt)
} # Finished in 3mins 50secs with 8 cores 02Feb23

bc.count.df <- as.data.frame(bind_rows(a))
#create bc.count.csv
bc.count.df[is.na(bc.count.df)] <- 0
bc.count.df <- bc.count.df %>% group_by(barcodes) %>% summarise_each(funs(sum))
bc.count.df <- bc.count.df %>% separate(col = barcodes, into = c("barcodes","hORF"),sep = "-")
# mipsa.metadata <- read_excel(path = "~/data-hlarman1/MIPSA-pep_db/RawData/MIPSApep_001/UCI.Pepseq.dictionary/20220919_MIPSA-Pep_AssDev_O-pool_annot.xlsx", col_names = TRUE)
# bc.count.df <- left_join(x = bc.count.df, y = mipsa.metadata[,c("Pepseq","Protein_name")], by = "Pepseq")
# bc.count.df <- bc.count.df %>% relocate(Protein_name, .after = Pepseq)

write.csv(bc.count.df, file = "mono.poly.UCI.v1.files/mipsa22.23.bc.count.mono+poly.v1.csv", row.names = FALSE)

#create orf.count.csv
orf.count.df <- bc.count.df[,-1] %>% group_by(hORF) %>% summarise_each(funs(sum))

write.csv(orf.count.df, file = "mono.poly.UCI.v1.files/mipsa22.23.orf.count.mono+poly.v1.csv", row.names = FALSE)

# sum(bc.count.df[grep("TCTGTCTCAGTGAGTTTCAC",bc.count.df$barcodes),"SapI_pur_franken_3_S9"])

