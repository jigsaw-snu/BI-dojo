# Current Working Directory : .../Project_CL-Treatment/Codes
# Raw Data Directory : .../Project_CL-Treatment/Resources/Raw
# Figures Directory : .../Project_CL-Treatment/Resources/Figures
# RData Dierectory : .../Project_CL-Treatment/Resources/RData




# 0. Global Variables (i.e. Path)
RSRC <- "../Resources"
RAW_DATA <- file.path(RSRC, "Raw")
FIGURES <- file.path(RSRC, "Figures")
RDATA <- file.path(RSRC, "RData")




# 1. Make Count Data
library(readxl)
library(tidyverse)  # dplyr, tibble, ggplot2

count_data <- readxl::read_xlsx(file.path(RAW_DATA, "SPF_con_vs_cl.xlsx"),
                                col_names = TRUE)

count_data <- count_data %>% 
    dplyr::select(c('Gene_Symbol', 
                    colnames(count_data)[grep(pattern = "Read_Count", x = colnames(count_data))]))

count_data <- count_data[, order(colnames(count_data))]

colnames(count_data) <- c("Gene_Symbol", 
                          paste0("SPF-CL-", as.character(seq(6))),
                          paste0("SPF-Saline-", as.character(seq(8))))

count_data <- count_data %>% tibble::column_to_rownames(var = "Gene_Symbol")

count_data <- count_data[apply(X = count_data,
                               MARGIN = 1,
                               FUN = function(x) {all(x != 0)})
                         , ]




# 2. Make Meta Data
meta_data <- data.frame(Treatment = sapply(X = colnames(count_data),
                                           FUN = function(x) {strsplit(x, '-')[[1]][2]}),
                        row.names = colnames(count_data))

meta_data$Treatment <- as.factor(meta_data$Treatment)
meta_data$Treatment <- stats::relevel(meta_data$Treatment, ref = "Saline")




# 3. Make DESeq Object
library(DESeq2)

dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_data,
                                      colData = meta_data,
                                      design = ~Treatment)

dds <- DESeq2::DESeq(dds)




# 4. Sample QC
# - PCA
# - Sample by Sample Heatmap
library(ggrepel)
library(ComplexHeatmap)
library(RColorBrewer)

rlog <- DESeq2::rlog(dds, blind = TRUE)
rlog <- SummarizedExperiment::assay(rlog)

pca <- stats::prcomp(t(rlog),
                     center = TRUE,
                     scale. = TRUE)  # row : Sample , column : feature (genes ; Principal Components)
percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), digits = 2)
sd_ratio <- sqrt(percent_var[2] / percent_var[1])

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2])
pca_df <- merge(pca_df, meta_data, by = "row.names") %>% 
    tibble::column_to_rownames("Row.names")

ggplot2::ggplot(data = pca_df, 
                aes(x = PC1, y = PC2)) +
    ggplot2::geom_point(aes(color = Treatment, size = 5)) + 
    ggplot2::xlab(paste0("PC1 : ", percent_var[1], '%')) +
    ggplot2::ylab(paste0("PC2 : ", percent_var[2], '%')) +
    ggplot2::coord_fixed(ratio = sd_ratio) + 
    ggrepel::geom_text_repel(aes(label = rownames(pca_df)))

#### Different PCA method --> Different result? ####
# 1. DESeq2::plotPCA(rlog, intgroup = colnames(meta_data), returnData = TRUE)
# 2. stats::prcomp(t(assay(rlog)), center = FALSE, scale. = FALSE)
# 3. stats::prcomp(t(assay(rlog)), center = TRUE, scale. = FALSE)
# 4. stats::prcomp(t(assay(rlog)), center = FALSE, scale. = TRUE)
# 5. stats::prcomp(t(assay(rlog)), center = TRUE, scale. = TRUE)
# --> What is difference between DESeq2::plotPCA() and stats::prcomp() ?
# --> Do we need center = TRUE, scale. = TRUE option in stats::prcomp() ?


corr <- stats::cor(rlog)  # cor(t(rlog)) : Gene by Gene Correlation

ComplexHeatmap::pheatmap(mat = corr,
                         top_annotation = meta_data,
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         color = rev(RColorBrewer::brewer.pal(11, "RdBu")))


# remove outlier samples and repeat above works
drops <- c("SPF-Saline-2", "SPF-Saline-3", "SPF-Saline-4", "SPF-Saline-8",
           "SPF-CL-4", "SPF-CL-5", "SPF-CL-6")
count_data <- count_data[, !(colnames(count_data) %in% drops)]
meta_data <- meta_data[!(rownames(meta_data) %in% drops), , drop = FALSE]




# 5. Differential Expression Analysis
# - Volcano Plot
# - Gene by Sample Heatmap
library(apeglm)
library(writexl)
library(EnhancedVolcano)

lfc_res <- DESeq2::lfcShrink(dds,
                             coef = DESeq2::resultsNames(dds)[2],
                             type = "apeglm")

sig_lfc_res <- lfc_res %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    tibble::rownames_to_column(var = "Gene_Symbol")

writexl::write_xlsx(sig_lfc_res, path = file.path(RDATA, "sig_lfc_res_FC15.xlsx"))

norm_cnt <- counts(dds, normalized = TRUE) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Gene_Symbol")

sig_norm_cnt <- norm_cnt %>%
    dplyr::filter(Gene_Symbol %in% sig_lfc_res$Gene_Symbol)


# volcano plot
EnhancedVolcano::EnhancedVolcano(toptable = lfc_res,
                                 lab = rownames(lfc_res),
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 pCutoff = 0.05,
                                 FCcutoff = 1.5,
                                 pointSize = 2,
                                 labSize = 4,
                                 axisLabSize = 10)

# manually drawn valcano plot using ggplot2



# gene by sample heatmap
ComplexHeatmap::pheatmap(mat = sig_norm_cnt %>% tibble::column_to_rownames(var = "Gene_Symbol"),
                         name = "Normalized Expression",
                         cluster_cols = TRUE,
                         cluster_rows = TRUE,
                         show_colnames = TRUE,
                         show_rownames = FALSE,
                         scale = "row",
                         color = rev(RColorBrewer::brewer.pal(11, "RdBu")))


# marker gene heatmap
# marker genes : https://www.frontiersin.org/files/Articles/599134/fendo-12-599134-HTML/image_m/fendo-12-599134-g001.jpg
brown_markers <- c("Slc27a1", "Nrg4", "P2rx5", "Slc36a2", "Tmem120b",
                   "Adrb3", "Cidea", "Clstn3", "Cox7a1", "Cox7a2", "Epsti1",
                   "Fgf21", "Hspb7", "Kcnk3", "Lhx8", "Mtus1", "Ppargc1a",
                   "Plin5", "Prdm16", "Prdm16os", "Sirt1", "Rps19bp1", "Ucp1",
                   "Mpzl2", "Eva1a", "Eva1b", "Eva1c","Eva1",
                   "Bmp7", "Ebf2", "Ednrb", "Mir133b", "Mir206", "Myf5", "Pdk4", 
                   "Prex1", "Zic1")
beige_markers <- c("Slc27a1", "Nrg4", "P2rx5", "Slc36a2", "Tmem120b",
                   "Adrb3", "Cidea", "Clstn3", "Cox7a1", "Cox7a2", "Epsti1",
                   "Fgf21", "Hspb7", "Kcnk3", "Lhx8", "Mtus1", "Ppargc1a",
                   "Plin5", "Prdm16", "Prdm16os", "Sirt1", "Rps19bp1", "Ucp1",
                   "Ebf3", "Mapre3", "Hoxc8", "Hoxc9", "Pdgfra", "Asc1", "Slc7a10",
                   "Aqp7", "Trip4", "Car4", "Tnfrsf9","Cd40", "Cited1", "Ear2", 
                   "Nr2f6", "Shox2", "Slc27a1", "Sp100", "Tbx1", "Litaf", "Tmem26")
white_markers <- c("Slc27a1", "Nrg4", "P2rx5", "Slc36a2", "Tmem120b",
                   "Ebf3", "Mapre3", "Hoxc8", "Hoxc9", "Pdgfra", "Asc1", "Slc7a10",
                   "Mpzl2", "Eva1a", "Eva1b", "Eva1c","Eva1",
                   "Fabp4", "Fbxo31", "Fbxo5", "Lep", "Lpl", "Nr1h3", "Nrip1", 
                   "Rb1", "Rbl1", "Retn", "Serpina3k", "Serpina3c", "Serpina3m",
                   "Tcf21", "Wfdc21")


browning_markers <- setdiff(union(brown_markers, beige_markers), white_markers)
fat_markers <- union(union(brown_markers, beige_markers), white_markers)

ComplexHeatmap::pheatmap(mat = sig_norm_cnt %>%
                             dplyr::filter(Gene_Symbol %in% browning_markers) %>%
                             tibble::column_to_rownames(var = "Gene_Symbol") %>%
                             as.matrix(),
                         border_gp = grid::gpar(col = "black", lty = 1, lwd=1),
                         legend = ComplexHeatmap::Legend(at = c(-2, -1, 0, 1, 2), title_position = "topcenter"),
                         cluster_cols = TRUE,
                         cluster_rows = TRUE,
                         show_colnames = TRUE,
                         show_rownames = TRUE,
                         scale = "row",
                         color = rev(RColorBrewer::brewer.pal(11, "RdBu")))

ComplexHeatmap::pheatmap(mat = sig_norm_cnt %>%
                             dplyr::filter(Gene_Symbol %in% fat_markers) %>%
                             tibble::column_to_rownames(var = "Gene_Symbol"),
                         name = "Normalized Expression",
                         cluster_cols = TRUE,
                         cluster_rows = TRUE,
                         show_colnames = TRUE,
                         show_rownames = TRUE,
                         scale = "row",
                         color = rev(RColorBrewer::brewer.pal(11, "RdBu")))

pheatmap::pheatmap(mat = sig_norm_cnt %>%
                             dplyr::filter(Gene_Symbol %in% fat_markers) %>%
                             tibble::column_to_rownames(var = "Gene_Symbol"),
                         cluster_cols = TRUE,
                         cluster_rows = TRUE,
                         show_colnames = TRUE,
                         show_rownames = TRUE,
                         scale = "row",
                         color = rev(RColorBrewer::brewer.pal(11, "RdBu")))


# marker gene dot plot
lfc_res_markers <- sig_lfc_res %>% dplyr::filter(Gene_Symbol %in% fat_markers)
lfc_res_markers <- lfc_res_markers[order(lfc_res_markers$log2FoldChange), ]
lfc_res_markers$Gene_Symbol <- factor(lfc_res_markers$Gene_Symbol,
                                      levels = )

ggplot2::ggplot(data = sig_lfc_res %>%
                    dplyr::filter(Gene_Symbol %in% fat_markers), 
                aes(x = log2FoldChange, y = Gene_Symbol)) +
    ggplot2::scale_y_discrete() +
    ggplot2::geom_point(size = 10, aes(color = -log(padj))) +
    ggplot2::scale_color_gradient(low = "#1E88E5", high = "#E53935") + 
    ggplot2::theme_bw() +
    guides(guide_legend(title = element_text("test")))



# 6. Gene ID Conversion
library(biomaRt)

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "mmusculus_gene_ensembl")

bm <- biomaRt::getBM(mart = mart,
                     filters = "external_gene_name",  # input ID type
                     attributes = c("external_gene_name", "ensembl_gene_id", "entrezgene_id"),  # output ID type
                     values = sig_lfc_res$Gene_Symbol)  # input gene list

# change column names
# same column name with original dataset is needed when merging two datasets
colnames(bm) <- c("Gene_Symbol", "ENSMUSG", "Entrez")

# - check NULL-mapping
# you should remove all the NULL-mapped entries
bm <- subset(bm,
             !is.na(Gene_Symbol) & trimws(Gene_Symbol) != '' &
             !is.na(ENSMUSG) & trimws(ENSMUSG) != '' &
             !is.na(Entrez) & trimws(Entrez) != '')

# merge converted ID with original dataset
sig_lfc_res_ens <- merge(sig_lfc_res, bm, by = "Gene_Symbol")  # will include duplicates




# 7. Functional Analysis
# - GO / KEGG / GSEA (HALLMARK_DB) / STRING
# - Barplot / Dotplot
# Tutorial : https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
# Tutorial : https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
# Tutorial : https://bioinformatics-core-shared-training.github.io/Bulk_RNAseq_Course_Nov22/Bulk_RNAseq_Course_Base/Markdowns/12_Gene_set_testing.html
library(clusterProfiler)
library(org.Mm.eg.db)

# Over-representation analysis


# Gene Set Enrichment Analysis


# Barplot and Dotplot using DAVID result
up_kegg <- read.table(file = file.path(RDATA, "DAVID_UP_FC25_KEGG.txt"),
                      sep = '\t',
                      header = TRUE)

up_kegg <- up_kegg %>% 
    dplyr::filter(FDR < 0.05) %>%
    dplyr::select(Term, Count, FDR)

# do not decreasingly order --> we want high count term as higher position in y-axis
# if you decreasingly order, then term with highest count will have lowest position in y-axis
up_kegg <- up_kegg[order(up_kegg$Count), ]
up_kegg$Term <- factor(up_kegg$Term, levels = rev(up_kegg$Term))

up_kegg$Term <- sapply(X = up_kegg$Term,
                       FUN = function(x) {strsplit(x, ':')[[1]][2]})

# Barplot
ggplot2::ggplot(data = up_kegg, 
                aes(x = Count, y = Term))+
    ggplot2::theme_bw() +
    ggplot2::geom_bar(stat = "identity", 
                      aes(fill = -log(FDR))) +
    ggplot2::scale_fill_gradient(low = "#1E88E5", high = "#E53935")

# Dotplot
ggplot2::ggplot(data = up_kegg, 
                aes(x = Count, y = Term)) +
    ggplot2::geom_point(size = 6, 
                        aes(color = -log(FDR))) +
    ggplot2::scale_color_gradient(low = "#1E88E5", high = "#E53935") +
    ggplot2::theme_bw(base_size = 15)

# Dotplot with highlight on specific term (application)
ggplot2::ggplot(data = up_kegg, 
                aes(x = Count, y = Term)) +
    ggplot2::scale_y_discrete() +  # you must specify that y-axis is discrete (factors) to ggplot 
    ggplot2::geom_rect(aes(xmin = -Inf, 
                           xmax = Inf, 
                           ymin = which(Term == "Thermogenesis") - 0.2, 
                           ymax = which(Term == "Thermogenesis") + 0.2),
                       fill = "red",
                       alpha = 0.03) +  # plot geom_rect() first to make this highlight box as background of main plot
    ggplot2::geom_point(size = 10, 
                        aes(color = -log(FDR))) +
    ggplot2::scale_color_gradient(low = "#1E88E5", high = "#E53935",  # legend color
                                  name = "-log(FDR)",  # legend title
                                  labels = c('5', '10', '15'), breaks = c(5, 10, 15)) +  # legend ticks
    ggplot2::theme_bw(base_size = 20) +  # use base_size to consistently increase text size
    ggplot2::theme(legend.title = ggplot2::element_text(size = 15),  # legend title size
                   # remember that levels(up_kegg$Term) and up_kegg$Term are opposite direction (reversed)
                   # so, c("black", "red", rep("black", 4)) will result not what you want
                   axis.text.y = ggplot2::element_text(face = c(rep("plain", 4), "bold", "plain"),  # axis label color and thickness
                                                       color = c(rep("black", 4), "red", "black")),
                   # margin (gap) betwwen axis and axis-title
                   axis.title.x = ggplot2::element_text(margin = margin(t = 20,  # top margin
                                                                        b = 0,  # bottom margin
                                                                        r = 0,  # right margin
                                                                        l = 0)),  # left margin
                   axis.title.y = ggplot2::element_text(margin = margin(t = 0,
                                                                        b = 0,
                                                                        r = 30,
                                                                        l = 0)))

# Dotplot with highlight on neighboring terms
ggplot2::ggplot(data = up_kegg,
                aes(x = Count, y = Term)) + 
    ggplot2::scale_y_discrete() +
    # Direction between levels(up_kegg$Term) and up_kegg$Term is opposite
    # Thermogenesis is 2nd term in up_kegg$Term, since up_kegg$Term is decreasingly ordered by counts
    # Thermogenesis is 5th term in levels(up_kegg$Term)
    # You should find 
    ggplot2::geom_rect(aes(xmin = -Inf,
                           xmax = Inf,
                           ymin = which(Term == "PPAR signaling pathway") - 0.2,
                           ymax = which(Term == "Thermogenesis") + 0.2,),
                       fill = "red",
                       alpha = 0.03) + 
    ggplot2::geom_point(size = 10,
                        aes(color = -log(FDR))) + 
    ggplot2::scale_color_gradient(low = "#1E88E5", high = "#E53935") +
    ggplot2::theme_bw()

# Dotplot with highlights on distant terms
ggplot2::ggplot(data = up_kegg,
                aes(x = Count, y = Term)) + 
    ggplot2::scale_y_discrete() +
    ggplot2::geom_rect(aes(xmin = -Inf,
                           xmax = Inf,
                           ymin = which(Term == "Ovarian steroidogenesis") - 0.2,
                           ymax = which(Term == "Ovarian steroidogenesis") + 0.2),
                       fill = "red",
                       alpha = 0.03) + 
    ggplot2::geom_rect(aes(xmin = -Inf,
                           xmax = Inf,
                           ymin = which(Term == "Thermogenesis") - 0.2,
                           ymax = which(Term == "Thermogenesis") + 0.2),
                       fill = "red",
                       alpha = 0.03) + 
    ggplot2::geom_point(size = 10,
                        aes(color = -log(FDR))) + 
    ggplot2::scale_color_gradient(low = "#1E88E5", high = "#E53935") + 
    ggplot2::theme_bw()




# 8. 
