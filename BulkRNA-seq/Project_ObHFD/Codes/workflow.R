RSRC <- "../Resources"
RAW_DATA <- file.path(RSRC, "Raw")
FIGURES <- file.path(RSRC, "Figures")
RDATA <- file.path(RSRC, "RData")




# 1. Prepare count_data and meta data

# Get (Raw) Dataset
library(GEOquery)
library(tidyverse)

GEOquery::getGEOSuppFiles("GSE167264", baseDir = "../Resources/Raw")

raw_count <- read.table(file.path(RAW_DATA, "GSE167264/GSE167264_Ob_Mus_RNA_seq_counts.txt.gz"),
                        skip = 1,
                        header = TRUE,
                        sep = '\t')

# select data from Epydidymal fat, Liver, Skeletal muscle for this analysis
SelwithMultiKeyword <- function(data, 
                                keywords,
                                MARGIN,
                                force_all = FALSE) {
    # MARGIN 
    # 0 : return logical selection
    # 1 : return row wise selected data
    # 2 : return column wise selected data
    # 3 : return selected data from vector
    MARGIN = as.character(MARGIN)
    
    switch(MARGIN,
           '0' = target <- data,
           '1' = target <- rownames(data),
           '2' = target <- colnames(data),
           '3' = target <- data)
    
    
    if (force_all) {
        sel <- Reduce('&', lapply(X = keywords, FUN = function(x) grepl(x, target)))
    } else {
        sel <- Reduce('|', lapply(X = keywords, FUN = function(x) grepl(x, target)))
    }
    
    
    switch(MARGIN,
           '0' = return(sel),
           '1' = return(data[sel, ]),
           '2' = return(data[, sel]),
           '3' = return(data[sel]))
}

count_data <- SelwithMultiKeyword(data = raw_count,
                                  keywords = c("Gene", "Ep", "Li", "Sk"),
                                  MARGIN = 2,
                                  force_all = FALSE)

# reset colnames
# Use interpretable keywords in name for future labels in plots
# Try use single '-' for the ease of making meta data --> strsplit(x, "-")
colnames(count_data)[grepl("Gene", colnames(count_data))] <- "Ens_Gene"  # gene name

colnames(count_data)[SelwithMultiKeyword(data = colnames(count_data),
                                         keywords = c("Ep", "WT", "ND"),
                                         MARGIN = 0, 
                                         force_all = TRUE)] <- paste0("Fat-Wild-NCD-", seq(3))

colnames(count_data)[SelwithMultiKeyword(data = colnames(count_data),
                                         keywords = c("Ep", "WT", "HFD"),
                                         MARGIN = 0, 
                                         force_all = TRUE)] <- paste0("Fat-Wild-HFD-", seq(3))

colnames(count_data)[SelwithMultiKeyword(data = colnames(count_data),
                                         keywords = c("Ep", "ob", "ND"),
                                         MARGIN = 0, 
                                         force_all = TRUE)] <- paste0("Fat-Obese-NCD-", seq(3))

colnames(count_data)[SelwithMultiKeyword(data = colnames(count_data),
                                         keywords = c("Ep", "ob", "HFD"),
                                         MARGIN = 0, 
                                         force_all = TRUE)] <- paste0("Fat-Obese-HFD-", seq(3))

colnames(count_data)[SelwithMultiKeyword(data = colnames(count_data),
                                         keywords = c("Li", "WT", "ND"),
                                         MARGIN = 0, 
                                         force_all = TRUE)] <- paste0("Liver-Wild-NCD-", seq(3))

colnames(count_data)[SelwithMultiKeyword(data = colnames(count_data),
                                         keywords = c("Li", "WT", "HFD"),
                                         MARGIN = 0, 
                                         force_all = TRUE)] <- paste0("Liver-Wild-HFD-", seq(3))

colnames(count_data)[SelwithMultiKeyword(data = colnames(count_data),
                                         keywords = c("Li", "ob", "ND"),
                                         MARGIN = 0, 
                                         force_all = TRUE)] <- paste0("Liver-Obese-NCD-", seq(3))

colnames(count_data)[SelwithMultiKeyword(data = colnames(count_data),
                                         keywords = c("Li", "ob", "HFD"),
                                         MARGIN = 0, 
                                         force_all = TRUE)] <- paste0("Liver-Obese-HFD-", seq(3))

colnames(count_data)[SelwithMultiKeyword(data = colnames(count_data),
                                         keywords = c("Sk", "WT", "ND"),
                                         MARGIN = 0, 
                                         force_all = TRUE)] <- paste0("Muscle-Wild-NCD-", seq(3))

colnames(count_data)[SelwithMultiKeyword(data = colnames(count_data),
                                         keywords = c("Sk", "WT", "HFD"),
                                         MARGIN = 0, 
                                         force_all = TRUE)] <- paste0("Muscle-Wild-HFD-", seq(3))

colnames(count_data)[SelwithMultiKeyword(data = colnames(count_data),
                                         keywords = c("Sk", "ob", "ND"),
                                         MARGIN = 0, 
                                         force_all = TRUE)] <- paste0("Muscle-Obese-NCD-", seq(3))

colnames(count_data)[SelwithMultiKeyword(data = colnames(count_data),
                                         keywords = c("Sk", "ob", "HFD"),
                                         MARGIN = 0, 
                                         force_all = TRUE)] <- paste0("Muscle-Obese-HFD-", seq(3))

count_data <- count_data %>% tibble::column_to_rownames(var = "Ens_Gene")

# sorting column names will ease your work making meta data
colnames(count_data) <- colnames(count_data)[order(colnames(count_data))] 

count_data <- count_data[apply(count_data, 1, function(x) all(x > 5)), ]


# making meta data
meta_data <- data.frame(Tissue = sapply(X = colnames(count_data),
                                        FUN = function(x) strsplit(x, "-")[[1]][1],
                                        USE.NAMES = FALSE),
                        Genotype = sapply(X = colnames(count_data),
                                          FUN = function(x) strsplit(x, "-")[[1]][2],
                                          USE.NAMES = FALSE),
                        Diet = sapply(X = colnames(count_data),
                                      FUN = function(x) strsplit(x, "-")[[1]][3],
                                      USE.NAMES = FALSE))

meta_data$Group <- paste(meta_data$Tissue, 
                         meta_data$Genotype, 
                         meta_data$Diet, 
                         sep = '-')

rownames(meta_data) <- colnames(count_data)

meta_data$Group <- factor(meta_data$Group)
meta_data$Tissue <- factor(meta_data$Tissue)
meta_data$Genotype <- factor(meta_data$Genotype, levels = c("Wild", "Obese"))
meta_data$Diet <- factor(meta_data$Diet, levels = c("NCD", "HFD"))




# 2. Make DESeq Object
library(DESeq2)
library(ggrepel)
library(PCAtools)

dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_data,
                                      colData = meta_data,
                                      design = ~Tissue * Genotype * Diet)

dds <- DESeq2::DESeq(dds)


    

# 3. Sample QC
# since each group has only 3 samples, we can't remove outliers
rlog <- DESeq2::rlog(dds, blind = TRUE)
rlog <- SummarizedExperiment::assay(rlog)

pca <- stats::prcomp(t(rlog), center = TRUE, scale. = TRUE)

percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), digits = 2)
sd_ratio <- sqrt(percent_var[2] / percent_var[1])

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2])
pca_df <- merge(pca_df, meta_data, by = "row.names") %>%
    tibble::column_to_rownames(var = "Row.names")

ggplot2::ggplot(data = pca_df, aes(x = PC1, y = PC2)) +
    ggplot2::geom_point(size = 5, aes(color = Diet)) +
    ggplot2::xlab(paste0("PC1 : ", percent_var[1], '%')) +
    ggplot2::ylab(paste0("PC2 : ", percent_var[2], '%')) +
    ggplot2::coord_fixed(ratio = sd_ratio)

ggplot2::ggplot(data = pca_df %>% dplyr::filter(Tissue == "Fat"), 
                aes(x = PC1, y = PC2)) +
    ggplot2::geom_point(size = 5, aes(color = Diet)) +
    ggplot2::xlab(paste0("PC1 : ", percent_var[1], '%')) +
    ggplot2::ylab(paste0("PC2 : ", percent_var[2], '%')) +
    ggplot2::coord_fixed(ratio = sd_ratio) +
    ggrepel::geom_text_repel(aes(label = rownames(pca_df %>% dplyr::filter(Tissue == "Fat"))))

ggplot2::ggplot(data = pca_df %>% dplyr::filter(Genotype == "Obese"), 
                aes(x = PC1, y = PC2)) +
    ggplot2::geom_point(size = 5, aes(color = Diet)) +
    ggplot2::xlab(paste0("PC1 : ", percent_var[1], '%')) +
    ggplot2::ylab(paste0("PC2 : ", percent_var[2], '%')) +
    ggplot2::coord_fixed(ratio = sd_ratio) +
    ggrepel::geom_text_repel(aes(label = rownames(pca_df %>% dplyr::filter(Genotype == "Obese"))))

ggplot2::ggplot(data = pca_df %>% dplyr::filter(Tissue == "Fat"), 
                aes(x = PC1, y = PC2)) +
    ggplot2::geom_point(size = 5, aes(color = Diet)) +
    ggplot2::xlab(paste0("PC1 : ", percent_var[1], '%')) +
    ggplot2::ylab(paste0("PC2 : ", percent_var[2], '%')) +
    ggplot2::coord_fixed(ratio = sd_ratio) +
    ggrepel::geom_text_repel(aes(label = rownames(pca_df %>% dplyr::filter(Tissue == "Fat"))))




# 4. 