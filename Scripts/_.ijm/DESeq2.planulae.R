# necessary libraries
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(ggforce)
library(ggVennDiagram)
library(NbClust)
library(ComplexHeatmap)
library(edgeR)
library(DEGreport)
library(fgsea)
library(purrr)
library(tibble)
library(reshape2)
library(patchwork)
library(stringr)
library(ggpattern)
library(inDAGO)
library(paletteer)


# setting working directory 
setwd("/home/gospozha/haifa/hiba/pl_align/")

#### Preparing necessary files ####

# list files with gene counts for each sample
dir = "/home/gospozha/haifa/hiba/pl_align/"
files = list.files(paste0(dir, "count.tables"), "*ReadsPerGene.out.tab", full.names = T)
files = list.files(paste0(dir, "count.tables.old"), "*ReadsPerGene.out.tab", full.names = T)
countData = data.frame(fread(files[1]))[c(1,3)]

# looping and reading the 2nd column from the remaining files
for(i in 2:length(files)) {
  countData = cbind(countData, data.frame(fread(files[i]))[3])
}

# skipping the first 4 lines, since count data starts on the 5th line
countData = countData[c(5:nrow(countData)),]

# renaming columns as sample names
#colnames(countData) = c("GeneID", gsub(paste0(dir,"count.tables/"), "", files))
# Get just the file names without the path
fn <- basename(files)

# Strip the suffix to get sample names (adjust pattern if needed)
sample_names <- sub("_ReadsPerGene\\.out\\.tab$", "", fn)

# Now set colnames
colnames(countData) <- c("GeneID", sample_names)
colnames(countData) = gsub("ReadsPerGene.out.tab", "", colnames(countData))
rownames(countData) = countData$GeneID
countData = countData[,c(2:ncol(countData))]
names <- colnames(countData)

# writing count matrix to a file
write.csv(countData, file="CountMatrix.csv")  
write.csv(countData, file="CountMatrix.old.csv")  
# reading count matrix from a file
countData  <- read.csv2('CountMatrix.csv', header=TRUE, row.names=1, sep=',', check.names = F)
countData  <- read.csv2('CountMatrix.old.csv', header=TRUE, row.names=1, sep=',', check.names = F)
# reading metadata file
MetaData <- read.csv2('Metadata', header=TRUE, sep="\t")

MetaData$condition <- as.factor(MetaData$condition)
#MetaData$origin <- as.factor(MetaData$sample)

MetaData <- MetaData[match(colnames(countData), MetaData$id), ]



#### Initial quality check ####

# Convert counts to DGEList
dge <- DGEList(counts = countData)
dim(dge)
# [1] 27083    10
# remove low counts
smallestGroupSize <- 3
keep <- rowSums(dge$counts >= 5) >= smallestGroupSize
dge <- dge[keep,]
dim(dge)
# [1] 7143   10
# 10286 with 5 threshold
# 9896 with stranded data
# Calculate FPM (Fragments Per Million)
fpm_values <- cpm(dge, normalized.lib.sizes = TRUE)  # edgeR's CPM is equivalent to FPM

# Convert to long format for plotting
fpm_df <- as.data.frame(fpm_values) %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "FPM") %>%
  left_join(MetaData, by = c("Sample" = "id"))  # Merge with metadata

ggplot(fpm_df, aes(x = FPM, color = condition)) +
  geom_density(alpha = 0.3) +
  scale_x_log10() +
  theme_minimal() +
  labs(title="Density Plot of FPM Values per Condition",
       x="FPM (log10 scaled)")

# statistical comparison
anova_res <- aov(FPM ~ condition, data = fpm_df)
summary(anova_res)
TukeyHSD(anova_res)
# they are the same

# reads per samples
library_sizes <- colSums(countData)
library_sizes

# Plot
barplot(library_sizes,
        las=2,
        main="Library sizes (tag-seq)",
        ylab="Total reads")
abline(h = 5e5, col="red", lty=2)  # warning threshold

# low-count genes percentage
low_count_fraction <- apply(countData, 2, function(x) mean(x < 5))
low_count_fraction

barplot(low_count_fraction,
        las=2,
        main="Fraction of low-count genes (<10)",
        ylab="Fraction")

abline(h = 0.7, col="red", lty=2)
abline(h = 0.8, col="darkred", lty=2)

# mean-variance check

gene_means <- rowMeans(countData)
gene_vars  <- apply(countData, 1, var)

qplot(log10(gene_means + 1), log10(gene_vars + 1),
      alpha = 0.3,
      main="Mean–variance distribution (Tag-seq)")

# saturation curve

saturation_plot <- inDAGO:::Saturation(
  matrix = countData, 
  method = "sampling", # Specify your counting method if required by the package version
  max_reads = 30000000,   # Specify max reads if needed, adjust as appropriate
  palette = "RColorBrewer::Paired"          # Choose a color palette
)

# Display the plot
print(saturation_plot)

# low count samples and library sizes
lib <- colSums(countData)

# fraction of genes <10
frac_low10 <- apply(countData, 2, function(x) mean(x < 10))

data.frame(
  sample = colnames(countData),
  library_size = lib,
  frac_low10 = round(frac_low10, 3)
)

# Quick plot
par(mfrow = c(1,2))
barplot(lib, las=2, main = "Library sizes")
barplot(frac_low10, las=2, main="Fraction of genes <10", ylim=c(0,1))
abline(h=0.7, col="red", lty=2)

#### DESeq2 model ####
# creating DESeq2 object 
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = MetaData,
                              design = ~ condition)

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 5) >= smallestGroupSize
dds <- dds[keep,]
dim(dds) # 7143 genes have left

# running a model
dds <- DESeq(dds)
res <- results(dds)

# Plotting histograms of p-values
hist(res$pvalue, breaks=50, col="skyblue", main="~ condition",
     xlab="p-value", xlim=c(0,1), ylim=c(0, max(table(cut(res$pvalue, breaks=50)))))

# saving a DESeq2 model to an R object
# saveRDS(dds, file = "dds_site_rrna_rin_condition.rds")
#dds <- readRDS(file = "../dds_site_rrna_rin_condition.rds")

#### PCA and sample distances using rlog ####

# estimating size factors to determine if it's better to use rlog
# to transform our data. rlog is more robust if size factors differ a lot
SF <- estimateSizeFactors(dds) 
print(sizeFactors(SF))

# the same using rlog transformation
rlog <- rlog(dds)
mat <- assay(rlog)
# PCA plot
pcaData <- plotPCA(rlog, intgroup=c("condition"), ntop = 500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#pdf("PCA.full.pdf",width=7)
pca<-ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  ggtitle("PCA of gene counts") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw()
pca
ggsave("pca.jpg", pca, width = 6.5, height = 6)
norm.counts <- assay(rlog)
write.csv(norm.counts, file="./rlog.counts.csv")

#dev.off()

# sample distances
sampleDists <- dist(t(assay(rlog)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rlog$condition)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#pdf("Dist.all.pdf",width=7)
dist <- pheatmap(sampleDistMatrix,
                 clustering_distance_rows=sampleDists,
                 clustering_distance_cols=sampleDists,
                 col=colors)
dist
ggsave("dist.jpg", dist, width = 6, height = 6)
#dev.off()


res <- results(dds, contrast=c("condition","High_fluo","Non_fluo"))
summary(res)

res.ordered <- data.frame(res) %>%
  filter(padj<.05 & abs(log2FoldChange)>1)  %>%
  arrange(padj) %>%
  mutate(Expression = case_when(log2FoldChange > log(1) ~ "Fluorescent",
                                log2FoldChange < -log(1) ~ "Non-fluorescent"))

View(res.ordered)


head(res.ordered)


anno<- read.csv2("sp_genes.tsv", sep="\t", header = T)
anno$NCBI.GeneID <- as.character(anno$NCBI.GeneID)

geneid <- sub("LOC", "", rownames(res.ordered))
res.ordered$NCBI.GeneID <- geneid

res_annot <- res.ordered %>%
  left_join(anno, by = "NCBI.GeneID")

#mutate(gene_id = sub("prefix_", "", gene_id))

write.csv(res_annot, "res_annot.csv")



#### Spis ####

head(res.ordered)


anno<- read.csv2("spis_tabulated_annots.csv", sep=",", header = T)
res.ordered$geneid <- rownames(res.ordered)

res_annot <- res.ordered %>%
  left_join(anno, by = "geneid")

#mutate(gene_id = sub("prefix_", "", gene_id))

write.csv(res_annot, "res_annot.old.csv")




gfp_terms <- "(gfp|green fluorescent|fluorescent protein|chromophore|gfp-like|cp?gfp)"

gfp <- anno %>%
  filter(grepl("GFP|fluorescent|chromoprotein",
               Description,
               ignore.case = TRUE))


res_full <-  data.frame(res)
geneid <- sub("LOC", "", rownames(res_full))
res_full$NCBI.GeneID <- geneid
res_full_annot <- res_full %>%
  left_join(anno, by = "NCBI.GeneID")
View(res_full_annot)



## plotting


gfp_res <- res_annot %>%
  filter(grepl("GFP|fluorescent|chromoprotein",
               Description,
               ignore.case = TRUE))

bio_genes <- c(gfp_res$Symbol) 
bio_genes <- bio_genes[bio_genes %in% rownames(rlog)]

mat <- assay(rlog)[bio_genes, ]
mat_z <- t(scale(t(mat)))  # Z-score by gene
anno_col <- data.frame(condition = colData(rlog)$condition)
rownames(anno_col) <- colnames(rlog)  

write.csv(mat_z, "z-normalized_expression_gfp.csv")

# heatmap of selected genes from rlog
hitmap <- pheatmap(mat_z,
         annotation_col = anno_col,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE)  # optional for clean plots

ggsave("heatmap.jpg", hitmap, width = 6, height = 6)
# are these genes significantly participate in depth change?

lrt.biomin <- res.ordered[bio_genes, ]
lrt.biomin <- na.omit(lrt.biomin)

# visualizing the boxplots for these genes
sig_bio_genes <- rownames(lrt.biomin) 

# Make sure they exist in rld
sig_bio_genes <- sig_bio_genes[sig_bio_genes %in% rownames(rlog)]

# Extract rlog expression matrix for those genes
expr_mat <- assay(rlog)[sig_bio_genes, ]

# Transpose and convert to data.frame
df <- as.data.frame(t(expr_mat))
df$sample <- rownames(df)
df$condition <- colData(rlog)$condition[match(df$sample, rownames(colData(rlog)))]

# Pivot longer for ggplot
df_long <- df %>%
  pivot_longer(cols = all_of(sig_bio_genes),
               names_to = "gene",
               values_to = "expression")
write.csv(df_long, "expression_boxplots.csv")

# Plot
boxplot <- ggplot(df_long, aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  facet_wrap(~ gene, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(title = "rlog Expression of GFP genes",
       y = "rlog Expression",
       x = "Condition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave("boxplot.jpg", boxplot, width = 6, height = 6)

  
