############################# Spatial transcriptomics - GeoMx DSP #############################################

###1. INSTALL PACKAGES AND LIBRARIES:

library(standR)
library(SpatialExperiment)
library(limma)
library(ExperimentHub)
library(tidyverse)
library(readr)
library(NanoStringNCTools)
library(GeomxTools)
library(ggpubr)
library(edgeR)
library(statmod)
library(ggrepel)
library(DT)
library(msigdb)
library(GSEABase)
library(SpatialDecon)
library(vissE)
library(igraph)
library(scater)
library(tidyverse)
library(gtable)
library(variancePartition)
library(dplyr)
library(readr)
library(ggplot2)
library(scran)

###2. LOAD SEQUENCING DATA:

sampleAnnoFile <- read_tsv("C:/Users/Simona/OneDrive - Imperial College London/Desktop/Segment.txt") %>% as.data.frame()
featureAnnoFile <- read_tsv("C:/Users/Simona/OneDrive - Imperial College London/Desktop/2.txt") %>% as.data.frame()
countFile <- read_tsv("C:/Users/Simona/OneDrive - Imperial College London/Desktop/Count.txt") %>% as.data.frame()
countFile <- subset(countFile, select = -GenomeBuild)


###3. DATA PROCESSING:
#3.1. Add more info to dataframe:

sampleAnnoFile$Batch <- c(rep("1", 56), rep("2", 39))
sampleAnnoFile$Group <- c(rep("C", 56), rep("A", 17), rep("B", 22))
sampleAnnoFile$Case <- c(rep("PDC83", 14), rep("PDC110", 14), rep("C088", 14), rep("PDC68", 14), rep("PDC05", 9), rep("C041", 8), rep("PDC84", 11), rep("PDC59", 11))
sampleAnnoFile$Age <- c(rep("84", 14), rep("94", 14), rep("94", 14), rep("95", 14), rep("58", 9), rep("54", 8), rep("83", 11), rep("82", 11))
sampleAnnoFile$Sex <- c(rep("Female", 14), rep("Female", 14), rep("Male", 14), rep("Female", 14), rep("Male", 9), rep("Male", 8), rep("Male", 11), rep("Female", 11))
sampleAnnoFile$ROI <- c(
  'Cortex tau and amyloid beta',
  'Cortex amyloid beta', 
  'Cortex tau',
  'Cortex tau',
  'Cortex tau and amyloid beta',
  'Cortex amyloid beta',
  'Cortex tau and amyloid beta',
  'Cortex tau and amyloid beta',
  'Cortex tau and amyloid beta',
  'Cortex tau and amyloid beta',
  'Cortex tau',
  'Cortex tau',
  'Cortex tau',
  'Cortex tau',
  'Cortex tau',
  'Cortex amyloid beta',
  'Cortex tau',
  'Cortex tau and amyloid beta',
  'Cortex tau and amyloid beta',
  'Cortex tau',
  'Cortex tau and amyloid beta',
  'Cortex tau',
  'Cortex tau',
  'Cortex tau and amyloid beta',
  'Cortex amyloid beta',
  'Cortex tau',
  'Cortex tau and amyloid beta',
  'Cortex tau and amyloid beta',
  'Cortex tau and amyloid beta',
  'Cortex tau and amyloid beta',
  'Cortex tau and amyloid beta',
  'Cortex tau and amyloid beta',
  'Cortex tau and amyloid beta',
  'Cortex tau',
  'Cortex tau',
  'Cortex tau',
  'Cortex amyloid beta',
  'Cortex amyloid beta',
  'Cortex tau',
  'Cortex tau',
  'Cortex tau',
  'Cortex tau and amyloid beta',
  'Cortex tau and amyloid beta',
  'Cortex tau',
  'Cortex tau',
  'Cortex tau and amyloid beta',
  'Cortex tau',
  'Cortex amyloid beta',
  'Cortex tau and amyloid beta',
  'Cortex tau',
  'Cortex tau',
  'Cortex amyloid beta',
  'Cortex tau and amyloid beta',
  'Cortex tau and amyloid beta',
  'Cortex tau and amyloid beta',
  'Cortex tau',
  'Cortex control',
  'Cortex control',
  'Cortex control',
  'Cortex amyloid beta',
  'Cortex control',
  'Cortex control',
  'Cortex control',
  'White matter control',
  'White matter control',
  'Cortex control',
  'Cortex control',
  'Cortex control',
  'Cortex control',
  'Cortex control',
  'Cortex control',
  'White matter control',
  'White matter control',
  'Cortex control',
  'Cortex control',
  'Cortex control',
  'Cortex control',
  'Cortex tau and amyloid beta',
  'Cortex control',
  'Cortex control',
  'White matter control',
  'White matter control',
  'Cortex tau and amyloid beta',
  'Cortex tau',
  'Cortex tau',
  'Cortex control',
  'Cortex control',
  'Cortex control',
  'Cortex control',
  'White matter control',
  'White matter control',
  'Cortex control',
  'Cortex control',
  'Cortex tau',
  'Cortex tau')

sampleAnnoFile$ROI_short <- c(
  'tau and amyloid beta',
  'amyloid beta', 
  'tau',
  'tau',
  'tau and amyloid beta',
  'amyloid beta',
  'tau and amyloid beta',
  'tau and amyloid beta',
  'tau and amyloid beta',
  'tau and amyloid beta',
  'tau',
  'tau',
  'tau',
  'tau',
  'tau',
  'amyloid beta',
  'tau',
  'tau and amyloid beta',
  'tau and amyloid beta',
  'tau',
  'tau and amyloid beta',
  'tau',
  'tau',
  'tau and amyloid beta',
  'amyloid beta',
  'tau',
  'tau and amyloid beta',
  'tau and amyloid beta',
  'tau and amyloid beta',
  'tau and amyloid beta',
  'tau and amyloid beta',
  'tau and amyloid beta',
  'tau and amyloid beta',
  'tau',
  'tau',
  'tau',
  'amyloid beta',
  'amyloid beta',
  'tau',
  'tau',
  'tau',
  'tau and amyloid beta',
  'tau and amyloid beta',
  'tau',
  'tau',
  'tau and amyloid beta',
  'tau',
  'amyloid beta',
  'tau and amyloid beta',
  'tau',
  'tau',
  'amyloid beta',
  'tau and amyloid beta',
  'tau and amyloid beta',
  'tau and amyloid beta',
  'tau',
  'control',
  'control',
  'control',
  'amyloid beta',
  'control',
  'control',
  'control',
  'control',
  'control',
  'control',
  'control',
  'control',
  'control',
  'control',
  'control',
  'control',
  'control',
  'control',
  'control',
  'control',
  'control',
  'tau and amyloid beta',
  'control',
  'control',
  'control',
  'control',
  'tau and amyloid beta',
  'tau',
  'tau',
  'control',
  'control',
  'control',
  'control',
  'control',
  'control',
  'control',
  'control',
  'tau',
  'tau')
sampleAnnoFile$ROI_short_Group <- paste0(sampleAnnoFile$ROI_short, " ", sampleAnnoFile$Group)
sampleAnnoFile$ROI_Group <- paste0(sampleAnnoFile$ROI, " ", sampleAnnoFile$Group)

#3.2. Create spatial experiment object:

spe <- readGeoMx(countFile, 
                 sampleAnnoFile, 
                 featureAnnoFile,
                 colnames.as.rownames = c("TargetName", "SegmentDisplayName", "TargetName"),
                 coord.colnames = c("ROICoordinateX", "ROICoordinateY"),
)
dim(spe)

#3.3. Delete ROIs where aligned reads are equal to 0:

spe <- spe[, spe$AlignedReads != 0]
dim(spe)
spe$AlignedReads

#comment: from 95 ROIs selected, 51 ROIs remained. 18676 genes total

###4. QUALITY CONTROL:
#4.1. Gene level quality control (min count 2 because count 1 is for negative probes):

spe <- addPerROIQC(spe, min_count = 2,  rm_genes = TRUE)
dim(spe)
metadata(spe) |> names()

plotGeneQC(
  spe,
  top_n = 9, 
  ordannots = "Group", 
  col = Group, 
  point_size = 2)

#4.2. ROI level quality control (library size against nuclei count):

plotROIQC(spe, x_threshold = 1, color = Group)
qc <- colData(spe)$AOINucleiCount > 1
table(qc) #no ROIs removed
spe <- spe[, qc]

plotRLExpr(spe)
plotRLExpr(spe, ordannots = "Group", assay = 2, color = Group)

#4.3. Dimension reduction PCA:

drawPCA(spe, assay = 2, color = Group)

set.seed(100)
spe <- scater::runPCA(spe)
pca_results <- reducedDim(spe, "PCA")
drawPCA(spe, precomputed = pca_results, col = Group)

plotScreePCA(spe, precomputed = pca_results)
plotPairPCA(spe, col = Group, precomputed = pca_results, n_dimension = 4)
plotPCAbiplot(spe, n_loadings = 10, precomputed = pca_results, col = Group)

###5. DATA NORMALISATION USING TMM:

spe_tmm <- geomxNorm(spe, method = "TMM")
plotRLExpr(spe_tmm, assay = 2, color = Group) + ggtitle("TMM") #to look how normalisation helps remove technical variation
set.seed(100)
spe_tmm <- scater::runPCA(spe_tmm)
pca_results_tmm <- reducedDim(spe_tmm, "PCA")
plotPairPCA(spe_tmm, precomputed = pca_results_tmm, color = Group)

###6. BATCH EFFECT CORRECTION:
#6.1. Correction method: Remove Unwanted Variation 4 (RUV4)

spe_tmm <- findNCGs(spe_tmm, batch_name = "SlideName", top_n = 300)
metadata(spe_tmm) |> names()

for(i in seq(5)){
  spe_ruv <- geomxBatchCorrection(spe_tmm, factors = "Case",
                                  NCGs = metadata(spe_tmm)$NCGs, k = i)
  
  print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = Group, title = paste0("k = ", i)))
  
} #Best k=1

spe_ruv <- geomxBatchCorrection(spe_tmm, factors = "Case", 
                                NCGs = metadata(spe_tmm)$NCGs, k = 1)

set.seed(100)
spe_ruv <- scater::runPCA(spe_ruv)
pca_results_ruv <- reducedDim(spe_ruv, "PCA")
plotPairPCA(spe_ruv, precomputed = pca_results_ruv, color = Group, title = "RUV4, k = 1", n_dimension = 4)

#6.2. Compare raw data and RUV4 coerrected:

spe_list <- list(spe_tmm, spe_ruv)
plotClusterEvalStats(spe_list = spe_list,
                     bio_feature_name = "Case",
                     batch_feature_name = "SlideName",
                     data_names = c("Raw","RUV4"))

plotRLExpr(spe_tmm, assay = 2, color = Group) + ggtitle("Raw")
plotRLExpr(spe_ruv, assay = 2, color = Group) + ggtitle("RUV4") #Use RUV4 for further analyses

###7. DIFFERENTIAL EXPRESSION ANALYSIS:

library(limma)

#7.1. Convert to DGE list object:

dge <- SE2DGEList(spe_ruv)

#7.2. Include group and RUV4 to design matrix:

design <- model.matrix(~0  + Group + ruv_W1, data = colData(spe_ruv))
colnames(design)

#7.3. Tidy factor names:

colnames(design) <- gsub("Group", "",colnames(design))
colnames(design) <- gsub(" ","_",colnames(design))
colnames(design)

#7.4. Group comparisons:

contrast_BA <- makeContrasts(B - A, levels = colnames(design))
contrast_CA <- makeContrasts(C - A, levels = colnames(design))
contrast_CB <- makeContrasts(C - B, levels = colnames(design))


#7.5. Filter out genes with low coverage in the dataset to allow a more accurate mean-variance relationship and reduce the number of statistical tests (default 10)

keep <- filterByExpr(dge, design, min.count = 10)
table(keep)
rownames(dge)[!keep]
dge_all <- dge[keep, ]

#Comment: 617 genes left

#7.6. BCV check:

dge_all <- estimateDisp(dge_all, design = design, robust = TRUE)
plotBCV(dge_all, legend.position = "topleft", ylim = c(0, 1.3))
bcv_df <- data.frame(
  'BCV' = sqrt(dge_all$tagwise.dispersion),
  'AveLogCPM' = dge_all$AveLogCPM,
  'gene_id' = rownames(dge_all)
)
highbcv <- bcv_df$BCV > 0.8
highbcv_df <- bcv_df[highbcv, ]
points(highbcv_df$AveLogCPM, highbcv_df$BCV, col = "red")
text(highbcv_df$AveLogCPM, highbcv_df$BCV, labels = highbcv_df$gene_id, pos = 4)

#7.7. Group comparisons:

v <- voom(dge_all, design, plot = TRUE) 
fit <- lmFit(v)
fit_contrast_BA <- contrasts.fit(fit, contrasts = contrast_BA)
efit_BA <- eBayes(fit_contrast_BA, robust = TRUE)
results_efit_BA <- decideTests(efit_BA, FDR = 0.05)
summary_efit_BA <- summary(results_efit_BA)
summary_efit_BA
#Comment: Intermediate vs Double-Negative: 258 down-regulated, 232 not sign, 127 up-regulated

fit_contrast_CA <- contrasts.fit(fit, contrasts = contrast_CA)
efit_CA <- eBayes(fit_contrast_CA, robust = TRUE)
results_efit_CA <- decideTests(efit_CA, FDR = 0.05)
summary_efit_CA <- summary(results_efit_CA)
summary_efit_CA
#Comment: Double-Positive vs Double-Negative: 71 down-regulated; 503 not sign, 43 up-regulated

fit_contrast_CB <- contrasts.fit(fit, contrasts = contrast_CB)
efit_CB <- eBayes(fit_contrast_CB, robust = TRUE)
results_efit_CB <- decideTests(efit_CB, FDR = 0.05)
summary_efit_CB <- summary(results_efit_CB)
summary_efit_CB
#Comment: Double-Positive vs Intermediate: 124 down-regulated, 373 not sign, 120 up-regulated

#7.8. Visualise DEGs between groups (up and down)

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)

get_counts_from_summary <- function(s, label) {
  data.frame(
    comparison = label,
    down = as.integer(s["Down", 1]),
    up   = as.integer(s["Up", 1]),
    stringsAsFactors = FALSE
  )
}

counts <- bind_rows(
  get_counts_from_summary(summary_efit_CA, "C - A"),
  get_counts_from_summary(summary_efit_BA, "B - A"),
  get_counts_from_summary(summary_efit_CB, "C - B")
)

deg_counts_long <- counts %>%
  pivot_longer(c(up, down), names_to = "direction", values_to = "n_DEGs") %>%
  mutate(
    direction = recode(direction, up = "Up-regulated", down = "Down-regulated"),
    comparison = factor(comparison, levels = c("C - A", "B - A", "C - B"))
  )

ggplot(deg_counts_long, aes(x = comparison, y = n_DEGs, fill = direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(
    aes(label = n_DEGs),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 4
  ) +
  theme_bw() +
  labs(x = NULL, y = "Number of DEGs (FDR < 0.05)", fill = NULL) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    text = element_text(size = 14)
  )

de_results_BA <- topTable(efit_BA, coef = 1, n = Inf)
de_results_CA <- topTable(efit_CA, coef = 1, n = Inf)
de_results_CB <- topTable(efit_CB, coef = 1, n = Inf)

de_genes_toptable_BA <- topTable(efit_BA, coef = 1, sort.by = "P", n = Inf, p.value = 0.05)
de_genes_toptable_CA <- topTable(efit_CA, coef = 1, sort.by = "P", n = Inf, p.value = 0.05)
de_genes_toptable_CB <- topTable(efit_CB, coef = 1, sort.by = "P", n = Inf, p.value = 0.05)

#Double-Positive versus Intermediate:

de_results_CB <- de_results_CB %>%
  filter(!is.na(adj.P.Val))
de_results_CB <- de_results_CB %>%
  mutate(DE = ifelse(logFC > 0 & adj.P.Val < 0.05, "UP",
                     ifelse(logFC < 0 & adj.P.Val < 0.05, "DOWN", "NOT DE")))
table(de_results_CB$DE)

# Top 20 UP
top_up <- de_results_CB %>%
  filter(logFC > 0, adj.P.Val < 0.05) %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 20)

# Top 20 DOWN
top_down <- de_results_CB %>%
  filter(logFC < 0, adj.P.Val < 0.05) %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 20)

top_genes <- bind_rows(top_up, top_down) %>%
  rownames_to_column("gene")

ggplot(de_results_CB, aes(AveExpr, logFC, col = DE)) + 
  
  geom_point(shape = 1, size = 1) +
  
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.5,
    point.padding = 0.3
  ) +
  
  theme_bw() +
  xlab("Average log-expression") +
  ylab("Log-fold-change") +
  ggtitle("C vs B") +
  
  scale_color_manual(values = c(
    "UP" = "blue",
    "NOT DE" = "gray",
    "DOWN" = "red"
  )) +
  
  theme(text = element_text(size = 15))

#Double-Positive vs Double-Negative

de_results_CA <- de_results_CA %>%
  filter(!is.na(adj.P.Val))

de_results_CA <- de_results_CA %>%
  mutate(DE = ifelse(logFC > 0 & adj.P.Val < 0.05, "UP",
                     ifelse(logFC < 0 & adj.P.Val < 0.05, "DOWN", "NOT DE")))

table(de_results_CA$DE)

# Top 20 UP
top_up <- de_results_CA %>%
  filter(logFC > 0, adj.P.Val < 0.05) %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 20)

# Top 20 DOWN
top_down <- de_results_CA %>%
  filter(logFC < 0, adj.P.Val < 0.05) %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 20)

top_genes <- bind_rows(top_up, top_down) %>%
  rownames_to_column("gene")

ggplot(de_results_CA, aes(AveExpr, logFC, col = DE)) + 
  
  geom_point(shape = 1, size = 1) +
  
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.5,
    point.padding = 0.3
  ) +
  
  theme_bw() +
  xlab("Average log-expression") +
  ylab("Log-fold-change") +
  ggtitle("C vs A") +
  
  scale_color_manual(values = c(
    "UP" = "blue",
    "NOT DE" = "gray",
    "DOWN" = "red"
  )) +
  
  theme(text = element_text(size = 15))

#Intermediate vs Double-Negative:

de_results_BA <- de_results_BA %>%
  filter(!is.na(adj.P.Val))
de_results_BA <- de_results_BA %>%
  mutate(DE = ifelse(logFC > 0 & adj.P.Val < 0.05, "UP",
                     ifelse(logFC < 0 & adj.P.Val < 0.05, "DOWN", "NOT DE")))

table(de_results_CB$DE)

# Top 20 UP
top_up <- de_results_BA %>%
  filter(logFC > 0, adj.P.Val < 0.05) %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 20)

# Top 20 DOWN
top_down <- de_results_BA %>%
  filter(logFC < 0, adj.P.Val < 0.05) %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 20)

top_genes <- bind_rows(top_up, top_down) %>%
  rownames_to_column("gene")

ggplot(de_results_BA, aes(AveExpr, logFC, col = DE)) + 
  
  geom_point(shape = 1, size = 1) +
  
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.5,
    point.padding = 0.3
  ) +
  
  theme_bw() +
  xlab("Average log-expression") +
  ylab("Log-fold-change") +
  ggtitle("B VS A") +
  
  scale_color_manual(values = c(
    "UP" = "blue",
    "NOT DE" = "gray",
    "DOWN" = "red"
  )) +
  
  theme(text = element_text(size = 15))

### 8. CELLULAR DECONVULATION FOR GEOMX

#8.1. Upload packackaged and re-upload datasets:
BiocManager::install("SpatialDecon")
BiocManager::install("SingleCellExperiment")
BiocManager::install("SummarizedExperiment")

library(SpatialDecon)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(dplyr)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

sampleAnnoFile <- read_tsv("C:/Users/Simona/OneDrive - Imperial College London/Desktop/Segment.txt") %>% as.data.frame()
featureAnnoFile <- read_tsv("C:/Users/Simona/OneDrive - Imperial College London/Desktop/2.txt") %>% as.data.frame()
countFile <- read_tsv("C:/Users/Simona/OneDrive - Imperial College London/Desktop/Count.txt") %>% as.data.frame()
countFile <- subset(countFile, select = -GenomeBuild)

#Remove duplicate row names
countFile <- countFile[!duplicated(countFile$TargetName), ]
featureAnnoFile <- featureAnnoFile[!duplicated(featureAnnoFile$TargetName), ]

#Set row names for count data and feature annotation data
rownames(countFile) <- countFile$TargetName
rownames(featureAnnoFile) <- featureAnnoFile$TargetName

#Ensure NegProbe-WTX is present and unique
if (sum(grepl("NegProbe-WTX", rownames(countFile))) > 1) {
  countFile <- countFile[!duplicated(rownames(countFile)), ]
}
if (sum(grepl("NegProbe-WTX", rownames(featureAnnoFile))) > 1) {
  featureAnnoFile <- featureAnnoFile[!duplicated(rownames(featureAnnoFile)), ]
}

spe <- readGeoMx(countFile, 
                 sampleAnnoFile, 
                 featureAnnoFile,
                 colnames.as.rownames = c("TargetName", "SegmentDisplayName", "TargetName"),
                 coord.colnames = c("ROICoordinateX", "ROICoordinateY"),
                 rmNegProbe = FALSE
)


#Delete ROIs where aligned reads are equal to 0
spe <- spe[, spe$AlignedReads != 0]
dim(spe)
spe$AlignedReads

#8.2. Gene level quality control (min count 2 because count 1 is for negative probes)

spe <- addPerROIQC(spe, min_count = 2,  rm_genes = TRUE)
dim(spe)
metadata(spe) |> names()

qc <- colData(spe)$AOINucleiCount > 1  # Adjust the threshold as needed
spe <- spe[, qc]

spe_tmm <- geomxNorm(spe, method = "TMM")
spd <- prepareSpatialDecon(spe_tmm)

human_brain <- download_profile_matrix(species = "Human",
                                       age_group = "Adult", 
                                       matrixname = "Brain_AllenBrainAtlas")
dim(human_brain)
heatmap(sweep(human_brain, 1, apply(human_brain, 1, max), "/"),
        labRow = NA, margins = c(10, 5), cexCol = 0.7)

res <- spatialdecon(norm = spd$normCount,
                    bg = spd$backGround,
                    X = human_brain,
                    align_genes = TRUE)
head(res)
str(res)

#8.3. Visualise cell deconvulation

BiocManager::install("DESeq2")
library(dplyr)
library(DESeq2) 

#Extract the proportion matrix
prop_matrix <- res$prop_of_all

#Convert matrix to data frame and reshape it
prop_df <- as.data.frame(prop_matrix)
prop_df <- rownames_to_column(prop_df, "CellType")
prop_long <- prop_df %>%
  pivot_longer(-CellType, names_to = "Sample", values_to = "Proportion")

#Plot the data using ggplot2
ggplot(prop_long, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = .7) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "Cell Type Proportions Across Samples", x = "Sample", y = "Proportion")

#8.4. Differential proportion analysis

# Extract the proportion matrix
prop_matrix <- res$prop_of_all

# Convert matrix to data frame and reshape it
prop_df <- as.data.frame(prop_matrix)
prop_df <- rownames_to_column(prop_df, "CellType")
prop_long <- prop_df %>%
  pivot_longer(-CellType, names_to = "Sample", values_to = "Proportion")

# Creating col_data
col_data <- data.frame(
  Sample = colnames(prop_matrix),
  Group = c(rep("C", 16), rep("A", 13), rep("B", 22))
)

# Merge proportion data with group information
prop_long <- merge(prop_long, col_data, by.x = "Sample", by.y = "Sample")

library(lme4)

# Fit a linear mixed-effects model to test for differential proportions
fit <- lmer(Proportion ~ Group + (1|CellType), data = prop_long)

# Summarize the results
summary(fit)

# Bar plot
ggplot(prop_long, aes(x = Group, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = 0.7) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "Cell Type Proportions Across Groups", x = "Group", y = "Proportion")

celltype_map <- c(
  # Excitatory neurons
  "L2.3.IT" = "Excitatory neuron",
  "L4.IT" = "Excitatory neuron",
  "L5.IT" = "Excitatory neuron",
  "L5.ET" = "Excitatory neuron",
  "L6.IT" = "Excitatory neuron",
  "L6.CT" = "Excitatory neuron",
  "L6b" = "Excitatory neuron",
  "L5.6.IT.Car3" = "Excitatory neuron",
  "L5.6.NP" = "Excitatory neuron",
  
  # Inhibitory neurons
  "PVALB" = "Inhibitory neuron",
  "SST" = "Inhibitory neuron",
  "VIP" = "Inhibitory neuron",
  "LAMP5" = "Inhibitory neuron",
  
  # Glia
  "Astrocyte" = "Astrocyte",
  "Microglia" = "Microglia",
  "Oligodendrocyte" = "Oligodendrocyte",
  "OPC" = "OPC"
)
# Convert matrix to data frame
prop_df <- as.data.frame(prop_matrix)
prop_df <- rownames_to_column(prop_df, "CellType")

# Collapse cell types
prop_df$CellType_major <- celltype_map[prop_df$CellType]

# Keep only the 6 major cell types
prop_df <- prop_df %>% 
  filter(!is.na(CellType_major))

# Reshape
prop_long <- prop_df %>%
  select(-CellType) %>%
  rename(CellType = CellType_major) %>%
  pivot_longer(-CellType, names_to = "Sample", values_to = "Proportion")

col_data <- data.frame(
  Sample = colnames(prop_matrix),
  Group = c(rep("C", 16), rep("A", 13), rep("B", 22))
)

prop_long <- merge(prop_long, col_data, by = "Sample")
library(lme4)

fit <- lmer(Proportion ~ Group + (1 | CellType), data = prop_long)
summary(fit)
ggplot(prop_long, aes(x = Group, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack",
           color = "black", width = 0.7) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(size = 14)
  ) +
  labs(
    title = "Major Brain Cell Type Proportions Across Groups",
    x = "Group",
    y = "Proportion"
  )

library(limma)
library(dplyr)
library(tidyr)
library(ggplot2)

# Helper: extract significant genes + direction from decideTests output
extract_deg_table <- function(decide_obj, label) {
  tt <- as.data.frame(decide_obj)
  colnames(tt) <- "call"  # -1, 0, 1
  tt$Genes <- rownames(tt)
  
  tt %>%
    filter(call != 0) %>%
    mutate(
      comparison = label,
      direction = ifelse(call == 1, "Up-regulated", "Down-regulated")
    ) %>%
    select(comparison, Genes, direction)
}

deg_BA <- extract_deg_table(results_efit_BA, "B - A")
deg_CA <- extract_deg_table(results_efit_CA, "C - A")
deg_CB <- extract_deg_table(results_efit_CB, "C - B")

deg_all <- bind_rows(deg_BA, deg_CA, deg_CB)

gene_celltype_df <- Supplementary_Table_5_Differentially_Expressed_Genes_Across_Comparisons %>%
  select(Genes, cell_type) %>%
  filter(!is.na(Genes), !is.na(cell_type)) %>%
  distinct()

# Collapse to 6 major classes using flexible pattern matching
gene_celltype_df <- gene_celltype_df %>%
  mutate(
    CellType_major = case_when(
      # Excitatory: layer/class labels (IT/ET/CT/NP/Car3/L6b etc) OR "excit"
      grepl("^L[0-9]", cell_type) | grepl("IT|ET|CT|NP|Car3|L6b", cell_type) ~ "Excitatory neuron",
      grepl("excit", cell_type, ignore.case = TRUE) ~ "Excitatory neuron",
      
      # Inhibitory: canonical interneuron subclasses OR "inhib"/"interneuron"
      grepl("PVALB|SST|VIP|LAMP5", cell_type) ~ "Inhibitory neuron",
      grepl("inhib|interneuron", cell_type, ignore.case = TRUE) ~ "Inhibitory neuron",
      
      # Glia
      grepl("astro", cell_type, ignore.case = TRUE) ~ "Astrocyte",
      grepl("micro", cell_type, ignore.case = TRUE) ~ "Microglia",
      grepl("opc", cell_type, ignore.case = TRUE) ~ "OPC",
      grepl("oligo", cell_type, ignore.case = TRUE) ~ "Oligodendrocyte",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(CellType_major)) %>%
  select(Genes, CellType_major) %>%
  distinct()

# Sanity check: you should see all 6 major types here
print(table(gene_celltype_df$CellType_major))

#8.5. Join GeoMx DEGs to cell types and compute composition per comparison:

deg_with_celltype <- deg_all %>%
  left_join(gene_celltype_df, by = "Genes")

# Check mapping success
cat("Mapped DEG rows:", sum(!is.na(deg_with_celltype$CellType_major)), "of", nrow(deg_with_celltype), "\n")

# If this is low, inspect unmapped genes (often gene ID mismatch)
# head(deg_with_celltype %>% filter(is.na(CellType_major)) %>% distinct(Genes), 30)

# Keep only mapped
deg_with_celltype <- deg_with_celltype %>% filter(!is.na(CellType_major))

## ---- A) DEG ENTRY COUNTS (gene can appear multiple times across cell types) ----
celltype_counts <- deg_with_celltype %>%
  count(comparison, direction, CellType_major, name = "n_DEG_entries")

#Proportions within each comparison & direction
celltype_props <- celltype_counts %>%
  group_by(comparison, direction) %>%
  mutate(prop = n_DEG_entries / sum(n_DEG_entries)) %>%
  ungroup() %>%
  mutate(comparison = factor(comparison, levels = c("C - A", "B - A", "C - B")))

#Plot
celltype_counts <- celltype_counts %>%
  mutate(
    comparison = factor(comparison, levels = c("C - A", "B - A", "C - B"))
  )

ggplot(celltype_counts, aes(x = comparison, y = n_DEG_entries, fill = CellType_major)) +
  geom_col(color = "black", width = 0.75) +
  facet_wrap(~direction) +
  scale_y_continuous(
    limits = c(0, 600),
    breaks = seq(0, 600, 100),
    expand = c(0, 0)
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 20, hjust = 1),
    text = element_text(size = 14)
  ) +
  labs(
    x = NULL,
    y = "Number of DEG entries (FDR < 0.05)",
    fill = "Cell type"
  )

###9. GSEA with KEGG pathways
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
install.packages("org.Hs.eg.db")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

#9.1. KEGG pathways:

entrez_ids_CA <- bitr(de_genes_toptable_CA$HUGOSymbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
de_genes_CA <- merge(de_genes_toptable_CA, entrez_ids_CA, by.x = "HUGOSymbol", by.y = "SYMBOL")
gene_list_CA <- de_genes_CA$logFC
names(gene_list_CA) <- de_genes_CA$ENTREZID
gene_list_CA <- sort(gene_list_CA, decreasing = TRUE)

summary(gene_list_CA)

gsea_result_CA <- gseKEGG(
  geneList = gene_list_CA,
  organism = "hsa",
  minGSSize = 10,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

if (nrow(gsea_result_CA@result) == 0) {
  message("No pathways enriched under the specific pvalueCutoff.")
} else 
  # Create term similarity matrix
  gsea_result_CA <- pairwise_termsim(gsea_result_CA)

head(gsea_result_CA@result)
dotplot(gsea_result_CA, showCategory = 15) + ggtitle("Enriched KEGG Pathways C vs A")
emapplot(gsea_result_CA) + ggtitle("Enrichment Map of KEGG Pathways C vs A")

#Intermediate vs Double-Negative:
entrez_ids_BA <- bitr(de_genes_toptable_BA$HUGOSymbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
de_genes_BA <- merge(de_genes_toptable_BA, entrez_ids_BA, by.x = "HUGOSymbol", by.y = "SYMBOL")
gene_list_BA <- de_genes_BA$logFC
names(gene_list_BA) <- de_genes_BA$ENTREZID
gene_list_BA <- sort(gene_list_BA, decreasing = TRUE)

summary(gene_list_BA)

gsea_result_BA <- gseKEGG(
  geneList = gene_list_BA,
  organism = "hsa",
  minGSSize = 10,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

if (nrow(gsea_result_BA@result) == 0) {
  message("No pathways enriched under the specific pvalueCutoff.")
} else 
  # Create term similarity matrix
  gsea_result_BA <- pairwise_termsim(gsea_result_AB)

head(gsea_result_BA@result)
dotplot(gsea_result_BA, showCategory = 15) + ggtitle("Enriched KEGG Pathways: Intermediate versus Double-Negative")
emapplot(gsea_result_BA) + ggtitle("Enrichment Map of KEGG Pathways: Intermediate versus Double-Negative")

#Double-Positive vs Intermediate: 
entrez_ids_CB <- bitr(de_genes_toptable_CB$HUGOSymbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
de_genes_CB <- merge(de_genes_toptable_CB, entrez_ids_CB, by.x = "HUGOSymbol", by.y = "SYMBOL")
gene_list_CB <- de_genes_CB$logFC
names(gene_list_CB) <- de_genes_CB$ENTREZID
gene_list_CB <- sort(gene_list_CB, decreasing = TRUE)

summary(gene_list_CB)

gsea_result_CB <- gseKEGG(
  geneList = gene_list_CB,
  organism = "hsa",
  minGSSize = 10,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

if (nrow(gsea_result_CB@result) == 0) {
  message("No pathways enriched under the specific pvalueCutoff.")
} else 
  # Create term similarity matrix
  gsea_result_CB <- pairwise_termsim(gsea_result_CB)

head(gsea_result_CB@result)
dotplot(gsea_result_CB, showCategory = 15) + ggtitle("Enriched KEGG Pathways: Double-Positive vs Intermediate")
emapplot(gsea_result_CB) + ggtitle("Enrichment Map of KEGG Pathways: Double-Positive vs Intermediate")

#Double-Positive vs Double-Negative: 
entrez_ids_CA <- bitr(de_genes_toptable_CA$HUGOSymbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
de_genes_CA <- merge(de_genes_toptable_CA, entrez_ids_CA, by.x = "HUGOSymbol", by.y = "SYMBOL")
gene_list_CA <- de_genes_CA$logFC
names(gene_list_CA) <- de_genes_CA$ENTREZID
gene_list_CA <- sort(gene_list_CA, decreasing = TRUE)

summary(gene_list_CA)

gsea_result_CA <- gseKEGG(
  geneList = gene_list_CA,
  organism = "hsa",
  minGSSize = 10,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

if (nrow(gsea_result_CA@result) == 0) {
  message("No pathways enriched under the specific pvalueCutoff.")
} else 
  # Create term similarity matrix
  gsea_result_CA <- pairwise_termsim(gsea_result_CA)

head(gsea_result_CA@result)
dotplot(gsea_result_CA, showCategory = 15) + ggtitle("Enriched KEGG Pathways: Double-Positive vs Double-Negative")
emapplot(gsea_result_CA) + ggtitle("Enrichment Map of KEGG Pathways: Double-Positive vs Double-Negative")

#9.2 GSEA with GO BP pathways:

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(fgsea)

#Double-Positive vs Intermediate
entrez_ids_CB <- bitr(de_genes_toptable_CB$HUGOSymbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
de_genes_CB <- merge(de_genes_toptable_CB, entrez_ids_CB, by.x = "HUGOSymbol", by.y = "SYMBOL")
gene_list_CB <- de_genes_CB$logFC
names(gene_list_CB) <- de_genes_CA$ENTREZID
gene_list_CB <- sort(gene_list_CB, decreasing = TRUE)

summary(gene_list_CB)

gsea_go_result_CB <- gseGO(
  geneList = gene_list_CB,
  OrgDb = org.Hs.eg.db,
  ont = "BP",  
  minGSSize = 10,
  pvalueCutoff = 0.1,  # More lenient p-value cutoff
  verbose = FALSE
)

# Check for results
if (nrow(gsea_go_result_CB@result) == 0) {
  message("No pathways enriched under the specific pvalueCutoff.")
} else 
  # Create term similarity matrix
  gsea_go_result_CB <- pairwise_termsim(gsea_go_result_CB)

head(gsea_go_result_CB@result)
dotplot(gsea_go_result_CB, showCategory = 15) + ggtitle("Top 15 Enriched GO BP Pathways Double-Positive versus Intermediate")

library(dplyr)
library(ggplot2)

gsea_top15 <- gsea_go_result_CB@result %>%
  filter(p.adjust < 0.1) %>%
  arrange(p.adjust) %>%
  slice_head(n = 15) %>%
  mutate(Direction = ifelse(NES > 0, "Up-regulated", "Down-regulated"))

ggplot(
  gsea_top15,
  aes(
    x = NES,
    y = reorder(Description, NES),
    color = Direction,
    size = -log10(p.adjust)
  )
) +
  geom_point() +
  scale_color_manual(values = c("Up-regulated" = "blue", "Down-regulated" = "red")) +
  theme_bw() +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    color = NULL,
    size = "-log10(adj. p-value)"
  )

#Double-Positive vs Double-Negative (no MF, CC wew found)
gsea_go_result_CA <- gseGO(
  geneList = gene_list_CA,
  OrgDb = org.Hs.eg.db,
  ont = "BP",  # Options: "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
  minGSSize = 10,
  pvalueCutoff = 0.1,  # More lenient p-value cutoff
  verbose = FALSE
)

if (nrow(gsea_go_result_CA@result) == 0) {
  message("No pathways enriched under the specific pvalueCutoff.")
} else 
  # Create term similarity matrix
  gsea_go_result_CA <- pairwise_termsim(gsea_go_result_CA)

head(gsea_go_result_CA@result)
dotplot(gsea_go_result_CA, showCategory = 15) + ggtitle("Top 15 Enriched GO BP Pathways Double-Positive versus Double-Negative")

gsea_top15 <- gsea_go_result_CA@result %>%
  filter(p.adjust < 0.1) %>%
  arrange(p.adjust) %>%
  slice_head(n = 15) %>%
  mutate(Direction = ifelse(NES > 0, "Up-regulated", "Down-regulated"))

ggplot(
  gsea_top15,
  aes(
    x = NES,
    y = reorder(Description, NES),
    color = Direction,
    size = -log10(p.adjust)
  )
) +
  geom_point() +
  scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue")) +
  theme_bw() +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    color = NULL,
    size = "-log10(adj. p-value)"
  )

#Intermediate vs Double-Negative:
gsea_go_result_BA <- gseGO(
  geneList = gene_list_BA,
  OrgDb = org.Hs.eg.db,
  ont = "BP",  # Options: "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
  minGSSize = 10,
  pvalueCutoff = 0.1,  # More lenient p-value cutoff
  verbose = FALSE
)

if (nrow(gsea_go_result_BA@result) == 0) {
  message("No pathways enriched under the specific pvalueCutoff.")
} else 
  # Create term similarity matrix
  gsea_go_result_BA <- pairwise_termsim(gsea_go_result_BA)

head(gsea_go_result_BA@result)
dotplot(gsea_go_result_BA, showCategory = 15) + ggtitle("Top 15 Enriched GO BB Pathways Intermediate versus Double-Negative")

gsea_top15 <- gsea_go_result_BA@result %>%
  filter(p.adjust < 0.1) %>%
  arrange(p.adjust) %>%
  slice_head(n = 15) %>%
  mutate(Direction = ifelse(NES > 0, "Up-regulated", "Down-regulated"))

ggplot(
  gsea_top15,
  aes(
    x = NES,
    y = reorder(Description, NES),
    color = Direction,
    size = -log10(p.adjust)
  )
) +
  geom_point() +
  scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue")) +
  theme_bw() +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    color = NULL,
    size = "-log10(adj. p-value)"
  )

###10. OVERLAP SNRNA-SEQ VS GEOMX

#10.1 unique DEGs in snRNA-seq:

deg_counts <- Supplementary_Table_5_Differentially_Expressed_Genes_Across_Comparisons %>%
  filter(!is.na(Genes), !is.na(fdr)) %>%
  filter(fdr < 0.05) %>%                      # set your threshold
  mutate(Direction.of.Effect = tolower(Direction.of.Effect)) %>%
  filter(Direction.of.Effect %in% c("up", "down")) %>%
  distinct(comparison, Direction.of.Effect, Genes) %>%  # <-- ensures uniqueness
  count(comparison, Direction.of.Effect, name = "n_genes") %>%
  complete(comparison, Direction.of.Effect = c("up", "down"),
           fill = list(n_genes = 0))

ggplot(deg_counts, aes(x = comparison, y = n_genes, fill = Direction.of.Effect)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_text(
    aes(label = n_genes),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 4
  ) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Number of unique DEGs (FDR < 0.05)",
    fill = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 14)
  )

#10.2 All DEG in snRNA-seq:

deg_counts <- Supplementary_Table_5_Differentially_Expressed_Genes_Across_Comparisons %>%
  filter(!is.na(Genes), !is.na(fdr)) %>%
  filter(fdr < 0.05) %>%                      
  mutate(Direction.of.Effect = tolower(Direction.of.Effect)) %>%
  filter(Direction.of.Effect %in% c("up", "down")) %>%
  # distinct(...)  <-- REMOVE THIS LINE
  count(comparison, Direction.of.Effect, name = "n_genes") %>%
  complete(
    comparison,
    Direction.of.Effect = c("up", "down"),
    fill = list(n_genes = 0)
  )

ggplot(deg_counts, aes(x = comparison, y = n_genes, fill = Direction.of.Effect)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_text(
    aes(label = n_genes),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 4
  ) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Number of DEGs (FDR < 0.05)",  # ← update label
    fill = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 14)
  )

#10.3. Unique genes in GeoMx

library(dplyr)
library(ggplot2)
library(VennDiagram)
library(UpSetR)

GeoMx_DEG_CB <- results_efit_CB[results_efit_CB[, 1] != 0, , drop = FALSE]
GeoMx_DEG_CB <- rownames(GeoMx_DEG_CB)

GeoMx_DEG_CA <- results_efit_CA[results_efit_CA[, 1] != 0, , drop = FALSE]
GeoMx_DEG_CA <- rownames(GeoMx_DEG_CA)

GeoMx_DEG_BA <- results_efit_BA[results_efit_BA[, 1] != 0, , drop = FALSE]
GeoMx_DEG_BA <- rownames(GeoMx_DEG_BA)


GeoMx_DEG_all <- unique(c(
  GeoMx_DEG_BA,
  GeoMx_DEG_CA,
  GeoMx_DEG_CB
))

#10.4 Unique genes in RNA-seq

DEG_genes_RNAseq <- Supplementary_Table_5_Differentially_Expressed_Genes_Across_Comparisons$Genes
DEG_genes_RNAseq <- unique(na.omit(DEG_genes_RNAseq))
anyDuplicated(DEG_genes_RNAseq)

#10.5 Overlap + uniques

common_genes <- intersect(GeoMx_DEG_all, DEG_genes_RNAseq)
common_genes_df <- as.data.frame(common_genes)
unique_to_GeoMx <- setdiff(GeoMx_DEG_all, DEG_genes_RNAseq)
unique_to_RNAseq <- setdiff(DEG_genes_RNAseq, GeoMx_DEG_all)

cat("Overlap:", length(common_genes), "\n")

venn.plot <- venn.diagram(
  x = list(
    "GeoMx" = GeoMx_DEG_all,
    "RNA-seq" = DEG_genes_RNAseq
  ),
  filename = NULL,
  fill = c("red", "blue"),
  alpha = 0.5,
  cat.col = c("red", "blue"),
  cat.pos = 0,
  cat.dist = 0.05,
  scaled = TRUE
)

#10.6 Plot the Venn diagram
grid.draw(venn.plot)

#Intermediate vs Double-Negative:
RNA_seq_BvsA <- Supplementary_Table_5_Differentially_Expressed_Genes_Across_Comparisons[Supplementary_Table_5_Differentially_Expressed_Genes_Across_Comparisons$comparison == "Intermediate_vs_Double-Negative", "Genes"]
print(RNA_seq_BvsA)
RNA_seq_BvsA_list <- RNA_seq_BvsA$Genes
common_genes_BvsA_RNA_seq <- intersect(GeoMx_DEG_BA, RNA_seq_BvsA$Genes)
common_genes_BvsA_RNA_seq

unique_to_your_deg_BvsA <- setdiff(GeoMx_DEG_BA, RNA_seq_BvsA_list)
unique_to_RNA_seq_deg_BvsA <- setdiff(RNA_seq_BvsA_list, GeoMx_DEG_BA)

cat("common_genes_BvsA_RNA_seq:\n")
print(common_genes_BvsA_RNA_seq) 
common_genes_BvsA_RNA_seq <- as.data.frame(common_genes_BvsA_RNA_seq)
venn.plot_BvsA <- venn.diagram(
  x = list(
    "GeoMx" = GeoMx_DEG_BA,
    "RNA-seq" = RNA_seq_BvsA_list
  ),
  filename = NULL,
  fill = c("red", "blue"),
  alpha = 0.5,
  cat.col = c("red", "blue"),
  cat.pos = 0,
  cat.dist = 0.05,
  scaled = TRUE
)

#Double-Positive vs Intermediate:
RNA_seq_CvsB <- Supplementary_Table_5_Differentially_Expressed_Genes_Across_Comparisons[Supplementary_Table_5_Differentially_Expressed_Genes_Across_Comparisons$comparison == "Double_Positive_vs_Intermediate", "Genes"]
print(RNA_seq_CvsB)
RNA_seq_CvsB_list <- RNA_seq_CvsB$Genes
common_genes_CvsB_RNA_seq <- intersect(GeoMx_DEG_CB, RNA_seq_CvsB$Genes)
common_genes_CvsB_RNA_seq

unique_to_your_deg_CvsB <- setdiff(GeoMx_DEG_CB, RNA_seq_CvsB_list)
unique_to_RNA_seq_deg_CvsB <- setdiff(RNA_seq_CvsB_list, GeoMx_DEG_CB)

cat("common_genes_CvsB_RNA_seq:\n")
print(common_genes_CvsB_RNA_seq) 
common_genes_CvsB_RNA_seq <- as.data.frame(common_genes_CvsB_RNA_seq)
venn.plot_CvsB <- venn.diagram(
  x = list(
    "GeoMx" = GeoMx_DEG_CB,
    "RNA-seq" = RNA_seq_CvsB_list
  ),
  filename = NULL,
  fill = c("red", "blue"),
  alpha = 0.5,
  cat.col = c("red", "blue"),
  cat.pos = 0,
  cat.dist = 0.05,
  scaled = TRUE
)

#Double-Positive vs Double-Negative:
RNA_seq_CvsA <- Supplementary_Table_5_Differentially_Expressed_Genes_Across_Comparisons[Supplementary_Table_5_Differentially_Expressed_Genes_Across_Comparisons$comparison == "Double_Positive_vs_Double_Negative", "Genes"]
print(RNA_seq_CvsA)
RNA_seq_CvsA_list <- RNA_seq_CvsA$Genes
common_genes_CvsA_RNA_seq <- intersect(GeoMx_DEG_CA, RNA_seq_CvsA$Genes)
common_genes_CvsA_RNA_seq

unique_to_your_deg_CvsA <- setdiff(GeoMx_DEG_CA, RNA_seq_CvsA_list)
unique_to_RNA_seq_deg_CvsA <- setdiff(RNA_seq_CvsA_list, GeoMx_DEG_CA)

cat("common_genes_CvsA_RNA_seq:\n")
print(common_genes_CvsA_RNA_seq) 
common_genes_CvsA_RNA_seq <- as.data.frame(common_genes_CvsA_RNA_seq)
venn.plot_CvsA <- venn.diagram(
  x = list(
    "GeoMx" = GeoMx_DEG_CA,
    "RNA-seq" = RNA_seq_CvsA_list
  ),
  filename = NULL,
  fill = c("red", "blue"),
  alpha = 0.5,
  cat.col = c("red", "blue"),
  cat.pos = 0,
  cat.dist = 0.05,
  scaled = TRUE
)

###11. COMPARE GENES WITH EXISTING DATASETS:

install.packages("VennDiagram")
install.packages("UpSetR")
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(UpSetR)

#11.1 DEGs from Kunkle et al (2019) dataset:

Kunkle <- c('CR1',
            'BIN1',
            
            'INPP5D', 
            
            'HLA-DRB1',
            
            'TREM2',
            
            'CD2AP',
            
            'NYAP1',
            
            'CLU', 
            
            'PTK2B',
            
            'EPHA1',
            
            'ECHDC3',
            
            'SP1',
            
            'MS4A2',
            
            'PICALM', 
            
            'SORL1',
            
            'SLC24A4', 
            
            'FERMT2',
            
            'ADAM10', 
            
            'IQCK', 
            
            'WWOX',
            
            'ACE',
            
            'ABCA7',
            
            'APOE', 
            
            'CASS4',
            
            'ADAMTS1')

Kunkle_df <- data.frame(Gene = Kunkle, stringsAsFactors = FALSE)

CvsB_list <- de_genes_toptable_CB$HUGOSymbol
CvsB_df <- data.frame(Gene = CvsB_list, stringsAsFactors = FALSE)

CvsA_list <- de_genes_toptable_CA$HUGOSymbol
CvsA_df <- data.frame(Gene = CvsA_list, stringsAsFactors = FALSE)

BvsA_list <- de_genes_toptable_BA$HUGOSymbol
BvsA_df <- data.frame(Gene = BvsA_list, stringsAsFactors = FALSE)

#11.2 Find common genes:

common_genes_all <- intersect(common_genes, Kunkle_df$Gene)
common_genes_cvsB <- intersect(common_genes_CvsB_RNA_seq, Kunkle_df$Gene)
common_genes_CvsA <- intersect(common_genes_CvsA_RNA_seq, Kunkle_df$Gene)
common_genes_BvsA <- intersect(common_genes_BvsA_RNA_seq$Gene, Kunkle_df$Gene)

cat("Common genes_all:\n")
print(common_genes_all) #CLU#
cat("Common genes_CvsB:\n")
print(common_genes_CvsB) #CLU#
cat("Common genes_CvsA:\n")
print(common_genes_CvsA) #NOTHING#
cat("Common genes_BvsA:\n")
print(common_genes_BvsA) #SORL1, BIN1, CLU#
