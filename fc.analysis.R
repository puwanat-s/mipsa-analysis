library(tidyverse)
library(limma)
library(edgeR)
library(gplots)
library(RColorBrewer)
library(heatmaply)
# https://bioconductor.github.io/BiocWorkshops/rna-seq-analysis-is-easy-as-1-2-3-with-limma-glimma-and-edger.html
# https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29 -> voom paper

#create normalized cpm counts from a bc.count file ----
count <- as.data.frame(read_csv(file = "mono.UCI.v1.files/mipsa22.23.bc.count.mono.v1.csv", col_names = TRUE))
count <- count %>% mutate(barcode.horf = paste(barcodes,hORF,sep = "-"))
count <- count %>% relocate(barcode.horf, .before = barcodes)
rownames(count) <- count[,1]
mipsa22.count.dge <- edgeR::DGEList(count[,3:14]) #create DGE object from the count matrix
mipsa23.count.dge <- edgeR::DGEList(count[,16:39])
mipsa22.count.dge.norm <- edgeR::calcNormFactors(mipsa22.count.dge, method = "TMM")
mipsa23.count.dge.norm <- edgeR::calcNormFactors(mipsa23.count.dge, method = "TMM")
# count.dge.norm.cpm <- edgeR::cpm(count.dge.norm, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 0.1)
mipsa22.count.dge.norm.cpm <- edgeR::cpm(mipsa22.count.dge.norm, normalized.lib.sizes = TRUE)
mipsa23.count.dge.norm.cpm <- edgeR::cpm(mipsa23.count.dge.norm, normalized.lib.sizes = TRUE)
mipsa22.keepers <- rowSums(mipsa22.count.dge.norm.cpm>5)>=3
mipsa23.keepers <- rowSums(mipsa23.count.dge.norm.cpm>3)>=2
mipsa22.count.dge.norm.cpm <- mipsa22.count.dge.norm.cpm[mipsa22.keepers,]
mipsa23.count.dge.norm.cpm <- mipsa23.count.dge.norm.cpm[mipsa23.keepers,]
#a <- edgeR::cpm(count.dge.norm, log = FALSE, prior.count = 0.5)

#read in study design
mipsa22.study_design <- read_tsv("MIPSA22_studydesign.txt")
mipsa23.study_design <- read_tsv("MIPSA23_studydesign.txt")
#study_design <- study_design[,1:2]
mipsa22.group <- factor(mipsa22.study_design$group)
mipsa23.group <- factor(mipsa23.study_design$group)
mipsa22.design <- model.matrix(~0+mipsa22.group)
mipsa23.design <- model.matrix(~0+mipsa23.group)

mipsa22.count.dge.norm.cpm <- mipsa22.count.dge.norm.cpm[,mipsa22.study_design$sample]
mipsa23.count.dge.norm.cpm <- mipsa23.count.dge.norm.cpm[,mipsa23.study_design$sample]
mipsa22.v.count.dge.norm.cpm <- voom(mipsa22.count.dge.norm.cpm , mipsa22.design, plot = TRUE)
mipsa23.v.count.dge.norm.cpm <- voom(mipsa23.count.dge.norm.cpm , mipsa23.design, plot = TRUE)


fit <- lmFit(mipsa23.v.count.dge.norm.cpm, mipsa23.design)
mipsa22.contrast.matrix <- limma::makeContrasts(Franken = mipsa22.group3 - mipsa22.group1,
                                                Franken_salmon = mipsa22.group4 - mipsa22.group2,
                                                levels=mipsa22.design) 
mipsa23.contrast.matrix <- limma::makeContrasts(V32_WT_G = mipsa23.group5 - mipsa23.group1,
                                                V32_WT_K = mipsa23.group6 - mipsa23.group2,
                                                V32_WT_S = mipsa23.group7 - mipsa23.group3,
                                                V32_WT_T = mipsa23.group8 - mipsa23.group4,
                                                V40_WT_G = mipsa23.group9 - mipsa23.group1,
                                                V40_WT_K = mipsa23.group10 - mipsa23.group2,
                                                V40_WT_S = mipsa23.group11 - mipsa23.group3,
                                                V40_WT_T = mipsa23.group12 - mipsa23.group4,
                                                V32_V40_G = mipsa23.group5 - mipsa23.group9,
                                                V32_V40_K = mipsa23.group6 - mipsa23.group10,
                                                V32_V40_S = mipsa23.group7 - mipsa23.group11,
                                                V32_V40_T = mipsa23.group8 - mipsa23.group12,
                                                levels=mipsa23.design)
fits <- contrasts.fit(fit, mipsa23.contrast.matrix)
ebFit <- eBayes(fits)
result <- limma::decideTests(ebFit, method = "global", adjust.method = "BH", p.value = 0.05, lfc = 1)
head(result)
summary(result)
tophits <- topTable(ebFit, coef = c(1:12), number = Inf, p.value = 0.05, lfc = 1)

count.dge.norm.cpm.log <- cpm(mipsa23.count.dge.norm, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 0.5)
count.dge.norm.cpm.log <- rownames_to_column(as.data.frame(count.dge.norm.cpm.log), var = "barcodes.horf")
diffGenes <- count.dge.norm.cpm.log[rownames(tophits),]
diffGenes <- rownames_to_column(as.data.frame(diffGenes), var = "barcodes-hORF")
# diffGenes <- count.dge.norm.cpm.log[result[,5] !=0, ]
# head(diffGenes)
# dim(diffGenes)

clustRows <- hclust(as.dist(1-cor(t(diffGenes), method = "pearson")), method = "complete")
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=1)

module.color <- rainbow(length(unique(module.assign)), start = 0.1, end=0.9)
module.color <- module.color[as.vector(module.assign)]

myheatcolors <- brewer.pal(name="RdBu", n=11)
heatmap.2(diffGenes,
          Rowv=as.dendrogram(clustRows),
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=rev(myheatcolors), scale = "row", labRow = NA,
          density.info = "none", trace = "none",
          cexRow=0.7, cexCol=1.2, margins=c(15,15), srtCol = 45)

dev.off()

#myheatcolors2 <- colorRampPalette(colors=c("blue","white","red"))(100)
#heatmaply(diffGenes,
          # colors = myheatcolors2,
          # Rowv=as.dendrogram(clustRows),
          # RowSideColors=module.color,
          # #showticklabels=c(FALSE,FALSE),
          # scale='row', col = rev(myheatcolors), 
          # file = "mipsapep002.bc.heatmaply.html")
dev.off()

names(module.color) <- names(module.assign)
module.assign.df <- as_tibble(as.list(module.assign))
module.assign.pivot <- pivot_longer(module.assign.df,
                                    cols = 1:30,
                                    names_to = "barcodes.horf",
                                    values_to = "module")
module.assign.pivot <- module.assign.pivot %>%
  mutate(moduleColor = case_when(
    module == 1 ~ "#FF0000",
    module == 2 ~ "#00FF00",
    module == 3 ~ "#0000FF",
    module == 4 ~ "#FFFF00",
    module == 5 ~ "#00FFFF",
    module == 6 ~ "#FF00FF",
    module == 7 ~ "#000000",
    module == 8 ~ "#FFFFFF",
    module == 9 ~ "#A52A2A",
    module == 10 ~ "#ADD8E6",
    module == 11 ~ "#FF69B4",
    module == 12 ~ "#7CFC00"))

dev.off()
ggplot(module.assign.pivot) +
  aes(module) +
  geom_bar(aes(fill=moduleColor)) +
  theme_bw()

module.pick <- 1
my.module <- tophits[names(module.assign[module.assign %in% module.pick]),1:12]
my.module <- my.module %>% rownames_to_column(var = "barcodes-horf")
my.module <- my.module %>% separate(col = `barcodes-horf`, into = c("barcodes","horf"),sep = "-")
write_csv(my.module,"mipsa22.bc.log2fc.module2.beads.beads+salmon.csv")

#create normalized cpm counts from a pep.count file ----

# count <- as.data.frame(read_csv(file = "pep.count.csv", col_names = TRUE))
# rownames(count) <- count[,1]
# count.dge <- edgeR::DGEList(count[,3:44]) #create DGE object from the count matrix
# count.dge.norm <- edgeR::calcNormFactors(count.dge, method = "TMM")
# count.dge.norm.cpm <- edgeR::cpm(count.dge.norm, prior.count = 0.5) #ignore log tranformation ! -> what does this do ??? 

count <- as.data.frame(read_csv(file = "mono.UCI.v1.files/mipsa22.23.orf.count.mono.v1.csv", col_names = TRUE))
rownames(count) <- count[,1]
mipsa22.count.dge <- edgeR::DGEList(count[,2:13]) #create DGE object from the count matrix
mipsa23.count.dge <- edgeR::DGEList(count[,14:37])

mipsa22.count.dge.norm <- edgeR::calcNormFactors(mipsa22.count.dge, method = "TMM")
mipsa23.count.dge.norm <- edgeR::calcNormFactors(mipsa23.count.dge, method = "TMM")
# count.dge.norm.cpm <- edgeR::cpm(count.dge.norm, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 0.1)
mipsa22.count.dge.norm.cpm <- edgeR::cpm(mipsa22.count.dge.norm, normalized.lib.sizes = TRUE)
mipsa23.count.dge.norm.cpm <- edgeR::cpm(mipsa23.count.dge.norm, normalized.lib.sizes = TRUE)
keepers1 <- rowSums(mipsa22.count.dge.norm.cpm>1)>=3
keepers2 <- rowSums(mipsa23.count.dge.norm.cpm>2)>=2
mipsa22.count.dge.norm.cpm <- mipsa22.count.dge.norm.cpm[keepers1,]
mipsa23.count.dge.norm.cpm <- mipsa23.count.dge.norm.cpm[keepers2,]

#keepers1.df <- as.data.frame(keepers1)
#keepers1.df <- rownames_to_column(as.data.frame(keepers1.df), var = "hORF")
#mipsa22.count.dge.norm.cpm.check <- rownames_to_column(as.data.frame(mipsa22.count.dge.norm.cpm), var = "hORF")

#read in study design
mipsa22.study_design <- read_tsv("MIPSA22_studydesign.txt")
mipsa23.study_design <- read_tsv("MIPSA23_studydesign.txt")
#study_design <- study_design[,1:2]
mipsa22.group <- factor(mipsa22.study_design$group)
mipsa23.group <- factor(mipsa23.study_design$group)
mipsa22.design <- model.matrix(~0+mipsa22.group)
mipsa23.design <- model.matrix(~0+mipsa23.group)

mipsa22.count.dge.norm.cpm <- mipsa22.count.dge.norm.cpm[,mipsa22.study_design$sample]
mipsa23.count.dge.norm.cpm <- mipsa23.count.dge.norm.cpm[,mipsa23.study_design$sample]
mipsa22.v.count.dge.norm.cpm <- voom(mipsa22.count.dge.norm.cpm , mipsa22.design, plot = TRUE)
mipsa23.v.count.dge.norm.cpm <- voom(mipsa23.count.dge.norm.cpm , mipsa23.design, plot = TRUE)

fit <- lmFit(mipsa23.v.count.dge.norm.cpm, mipsa23.design)
mipsa22.contrast.matrix <- limma::makeContrasts(Franken = mipsa22.group3 - mipsa22.group1,
                                                Franken_salmon = mipsa22.group4 - mipsa22.group2,
                                        levels=mipsa22.design) 
mipsa23.contrast.matrix <- limma::makeContrasts(V32_WT_G = mipsa23.group5 - mipsa23.group1,
                                                V32_WT_K = mipsa23.group6 - mipsa23.group2,
                                                V32_WT_S = mipsa23.group7 - mipsa23.group3,
                                                V32_WT_T = mipsa23.group8 - mipsa23.group4,
                                                V40_WT_G = mipsa23.group9 - mipsa23.group1,
                                                V40_WT_K = mipsa23.group10 - mipsa23.group2,
                                                V40_WT_S = mipsa23.group11 - mipsa23.group3,
                                                V40_WT_T = mipsa23.group12 - mipsa23.group4,
                                                V32_V40_G = mipsa23.group5 - mipsa23.group9,
                                                V32_V40_K = mipsa23.group6 - mipsa23.group10,
                                                V32_V40_S = mipsa23.group7 - mipsa23.group11,
                                                V32_V40_T = mipsa23.group8 - mipsa23.group12,
                                                levels=mipsa23.design) 

fits <- contrasts.fit(fit, mipsa23.contrast.matrix)
ebFit <- eBayes(fits)
result <- limma::decideTests(ebFit, method = "global", adjust.method = "BH", p.value = 0.05, lfc = 2) # for mipsa22
# result <- limma::decideTests(ebFit, method = "global")
head(result)
summary(result)
tophits <- topTable(ebFit, coef = c(1:12), number = Inf, p.value = 0.05, lfc = 2)
tophits <- tophits %>% rownames_to_column(var = "hORF")

count.dge.norm.cpm.log <- cpm(mipsa23.count.dge.norm, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 5)
# count.dge.norm.cpm.log.check <- rownames_to_column(as.data.frame(count.dge.norm.cpm.log), var = "Pepseq")
diffGenes <- count.dge.norm.cpm.log[rownames(tophits),] # for mipsa22
diffGenes <- count.dge.norm.cpm.log
# https://support.bioconductor.org/p/100731/ suggests to use edgeR cpm rather than voom $E because of the fixed prior.count = 0.5 ? 
# diffGenes <- v.count.dge.norm.cpm$E[rownames(tophits),]
# head(diffGenes)
# dim(diffGenes)

clustRows <- hclust(as.dist(1-cor(t(diffGenes), method = "pearson")), method = "complete")
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)

module.color <- rainbow(length(unique(module.assign)), start = 0.1, end=0.9)
module.color <- module.color[as.vector(module.assign)]

myheatcolors <- brewer.pal(name="RdBu", n=11)
heatmap.2(diffGenes,
          Rowv=as.dendrogram(clustRows),
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=rev(myheatcolors), scale = "row", labRow = NA,
          density.info = "none", trace = "none",
          cexRow=1, cexCol=1.5, margins=c(15,25), srtCol = 45)
#dev.off()

names(module.color) <- names(module.assign)
module.assign.df <- as_tibble(as.list(module.assign))
module.assign.pivot <- pivot_longer(module.assign.df,
                                    cols = 1:1248,
                                    names_to = "hORF",
                                    values_to = "module")
module.assign.pivot <- module.assign.pivot %>%
  mutate(moduleColor = case_when(
    module == 1 ~ "#FF0099",
    module == 2 ~ "#FF9900"))

dev.off()
ggplot(module.assign.pivot) +
  aes(module) +
  geom_bar(aes(fill=moduleColor)) +
  theme_bw()


module.pick <- 2
my.module <- tophits[names(module.assign[module.assign %in% module.pick]),1:2]
my.module <- my.module %>% rownames_to_column(var = "hORF")
#my.module <- left_join(x = my.module, y = bc.count.df[,c("Pepseq","Protein_name")], by = "Pepseq")
#my.module <- my.module %>% relocate(Protein_name, .after = Pepseq)
write_csv(my.module,"mono.poly.UCI.v1.files/mipsa22.mono+polyUCI.v1.horf.log2fc.module2.beads.beads+salmon.csv")

