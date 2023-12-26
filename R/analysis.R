library(tidyverse)
library(tidyquant)
library(ggpubr)
library(pheatmap)
library(ComplexHeatmap)
library(biomaRt)
require(clusterProfiler)
require(org.Mm.eg.db)
keytypes(org.Mm.eg.db)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(plyranges)
library(DESeq2)

counts <- read_delim("data/counts.tabular", 
                     delim="\t",
                     col_types = c("c",rep("d",50))
                     )
counts <- counts[-c(56885:56888),]
counts <- as.data.frame(counts)

colnames(counts) <- gsub("filtered_|_S[0-9]+_L001_R1_001.fastq.gz_2","",colnames(counts))
colnames(counts)[1] <- "gene_id"

## load the gtf file used for mapping
## from url: https://support.bioconductor.org/p/112469/
gr <- read_gff(file.path("data","gencode.vM33.annotation.gtf.gz")) %>% select(gene_id, gene_name, transcript_id) %>% 
  as.data.frame()

gr <- gr %>% 
  filter(gene_id %in% counts$gene_id) %>% 
  distinct(gene_id, gene_name)
rownames(gr) <- gr$gene_id

## set rownames(counts)
rownames(counts) <- make.unique(gr[counts$gene_id,]$gene_name)
counts <- counts[,-1]
counts <- counts[,grepl("Lib_1130", colnames(counts))]
counts <- counts[,-8]
coldata <- data.frame(condition=case_when(grepl("Slice[1-9]_ctrl", colnames(counts))~"Slice_Ctrl",
                                          grepl("Ctrl", colnames(counts))~"Ctrl",
                                          grepl("IOMM", colnames(counts))~"IOMM-Lee",
                                          grepl("Myelin", colnames(counts))~"Myelin",
                                          T~"Slice"
                                          ),
                      batch=case_when(grepl("Lib_1130", colnames(counts))~ "B",
                                      grepl("Lib_1201", colnames(counts))~ "C",
                                      T~"A"
                                      ))
  
## run DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design= ~ condition)
## filter
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

## re-order
dds$condition <- factor(dds$condition, levels = c("Ctrl","Slice_Ctrl","Slice","IOMM-Lee","Myelin"))
dds$condition <- relevel(dds$condition, ref = "Ctrl")

## visualize and normalize
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

#normalized counts_genes
normalized_counts_genes <- counts(dds, normalized=TRUE)

#log transformation
vsd_wt <- vst(dds, blind=TRUE)

# Extract the vst matrix from the object
vsd_cor_wt <- vsd_wt %>% 
  assay() %>%
  cor()

# plot heatmap
pheatmap(vsd_cor_wt, annotation = select(coldata, condition),cellheight=8, cellwidth = 8)

pdf(file.path("plots","correlation_heatmap.pdf"))
pheatmap(vsd_cor_wt, annotation = select(coldata, condition),cellheight=45, cellwidth = 45)
dev.off()

# Plot PCA
plotPCA(vsd_wt, intgroup="condition") +
  scale_color_tq() +
  theme_linedraw() + 
  geom_point(size=7)

ggsave(file.path("plots","sample_pca.pdf"))

## run DESeq
dds <- DESeq(dds)

## plot dispersion
plotDispEsts(dds)

res <- results(dds)
res
res <- results(dds, name="condition_BCG_vs_CTRL")
res <- results(dds, contrast=c("condition","BCG","CTRL"))
resultsNames(dds)

## Log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef="condition_BCG_vs_CTRL", type="apeglm")
resLFC

plotMA(resLFC, ylim=c(-2,2))

idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

## significant results 
res05 <- results(dds, alpha=0.05)
summary(res05)

res_all <- data.frame(res05) %>%
  rownames_to_column(var = "Gene") 

#extract significant genes 
res_sig <- subset(res_all, padj < 0.05)
res_sig <- res_sig %>%
  arrange(log2FoldChange)

##write.csv()
#plot sig genes
# Subset normalized counts_genes to significant genes
sig_norm_counts_genes <- normalized_counts_genes[res_sig$Gene, ]

column_ha = HeatmapAnnotation(condition = coldata$condition)
ann_colors = list(
  condition = c(CTRL=unname(palette_light())[1], BCG=unname(palette_light())[2]))

# Run pheatmap
pheat <- pheatmap(sig_norm_counts_genes,
                  color = viridis(100, option = "C"),
                  cluster_rows = T,
                  show_rownames = T,
                  annotation_col =coldata,
                  annotation_colors = ann_colors,
                  scale = "row")
pheat

pdf(file.path("plots","diff_gene_heatmap.pdf"))
pheat
dev.off()

## define diffgenes table
diffgenes <- data.frame("gene" = res_sig$Gene,
                        "cluster" = case_when(
                          res_sig$log2FoldChange>0 ~ "BCG",
                          T ~ "CTRL"
                        ))
## plot
genes <- bitr(unique(diffgenes$gene), fromType = "SYMBOL",
              toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
              OrgDb = 'org.Mm.eg.db')

colnames(genes)[1] <- "gene"

diffgenes <- diffgenes %>% 
  left_join(genes, relationship = "many-to-many") %>% 
  distinct(cluster, gene, .keep_all = T) %>%
  group_by(cluster) %>% 
  #top_n(100, wt=avg_log2FC) %>% 
  na.omit()

diffgenes$cluster <- factor(diffgenes$cluster, levels = c("CTRL","BCG"))

go_terms_bp <- compareCluster(ENTREZID ~ cluster , 
                              data=diffgenes, 
                              fun = "enrichGO",
                              OrgDb='org.Mm.eg.db',
                              ont="BP"
)
dotplot(go_terms_bp)  
ggsave(file.path("plots","GO_term_BP_dotplot.pdf"), height = 7, width = 5)
