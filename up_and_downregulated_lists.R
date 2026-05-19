library(tidyverse)
library(tximport)
library(edgeR)
library(matrixStats)
library(cowplot)
library(Biostrings)
library(limma)
library(DESeq2)
library(dplyr)
library(tibble)

# set working directory to the folder where the study design file will be read
# and where output CSVs will be written
setwd("C:/Users/lisam/Documents/RNA_Seq_2026/gene_expression")

# path to the kallisto quantification folders
# each sample should have its own folder containing abundance.h5 or abundance.tsv
kallisto_dir <- "C:/Users/lisam/Documents/RNA_Seq_2026/kallisto_outputs_t2t_48hr_gapmer"

# fasta file used to build the kallisto index
# we use the fasta headers here to build a transcript-to-gene mapping
fasta_file <- "C:/Users/lisam/Documents/GapmerData/Lisa_M_Gapmer_2025/kallisto_t2t/t2t_censat.fna"

# read all fasta sequences so we can access the header lines
fa      <- readDNAStringSet(fasta_file)
headers <- names(fa)

# build a transcript to gene table from the fasta headers
# raw_id / tx_id: transcript identifier used by kallisto
# gene_id_satellite: tries to extract a satellite-style identifier after "censat_"
# gene_id_parens: fallback identifier from text inside parentheses
# gene_id: final chosen gene name used for gene-level DE testing
tx2gene_full <- tibble(header = headers) %>%
  mutate(
    raw_id = str_extract(header, "^\\S+"),
    tx_id  = raw_id,
    gene_id_satellite = str_extract(header, "(?<=censat_)[^ ]+"),
    gene_id_parens    = str_extract(header, "\\(([^)]+)\\)") %>% str_remove_all("[()]"),
    gene_id           = coalesce(gene_id_satellite, gene_id_parens, tx_id)
  ) %>%
  filter(!is.na(tx_id)) %>%
  distinct()

# keep only the columns DESeq2 / tximport need for gene-level summarization:
# transcript ID -> gene ID
tx2gene <- tx2gene_full %>% select(tx_id, gene_id)

# read metadata describing each sample
# treatment is converted to uppercase and stored as a factor so MOCK is the reference level
targets <- read_tsv("studydesign.txt", col_types = cols()) %>%
  mutate(
    condition = toupper(Treatment),
    replicate = as.character(Replicates)
  ) %>%
  mutate(
    condition = factor(condition, levels = c("MOCK", "M1", "SAT2"))
  )

# sample_names must match the kallisto output folder names
sample_names <- targets$Name

# create full paths to the kallisto abundance files
# prefer abundance.h5 if present, otherwise use abundance.tsv
h5_paths  <- file.path(kallisto_dir, sample_names, "abundance.h5")
tsv_paths <- file.path(kallisto_dir, sample_names, "abundance.tsv")
files <- ifelse(file.exists(h5_paths), h5_paths, tsv_paths)
names(files) <- sample_names

# stop immediately if any expected quantification files are missing
missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  stop("Missing abundance files:\n", paste(missing_files, collapse = "\n"))
}

# import kallisto estimates with tximport
# tx2gene collapses transcript-level estimates to gene-level counts
# countsFromAbundance = "lengthScaledTPM" adjusts counts using estimated abundance and transcript length
txi <- tximport(
  files,
  type = "kallisto",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM"
)

# prepare sample metadata for DESeq2
# row names of colData must match the column names of txi counts
colData <- targets %>%
  dplyr::select(Name, condition) %>%
  column_to_rownames("Name")

# create DESeq2 dataset
# design = ~ condition means gene expression is modeled as a function of treatment condition
dds <- DESeqDataSetFromTximport(
  txi,
  colData = colData,
  design = ~ condition
)

# explicitly set MOCK as the reference condition
# then log2 fold changes for SAT2 and M1 are interpreted relative to MOCK
dds$condition <- relevel(dds$condition, ref = "MOCK")

# run the full DESeq2 pipeline:
# 1. estimate size factors to normalize for library size
# 2. estimate gene-wise dispersion values
# 3. fit a negative binomial generalized linear model for each gene
# 4. perform Wald tests for the coefficients in the model
dds <- DESeq(dds)

# get normalized counts for descriptive summaries only
# these are not what DESeq2 directly tests, but are useful for reporting group averages
normc <- counts(dds, normalized = TRUE)


# calculate mean normalized count for each gene within each condition
# these values help interpret direction and magnitude of expression differences
group_means <- tibble(
  gene_name = rownames(normc),
  MOCK_avg  = rowMeans(normc[, dds$condition == "MOCK", drop = FALSE]),
  SAT2_avg  = rowMeans(normc[, dds$condition == "SAT2", drop = FALSE]),
  M1_avg    = rowMeans(normc[, dds$condition == "M1", drop = FALSE])
)

# -----------------------------
# SAT2 vs MOCK
# -----------------------------

res_sat2 <- results(
  dds,
  contrast             = c("condition", "SAT2", "MOCK"),
  alpha                = 0.05,
  pAdjustMethod        = "BH",
  independentFiltering = TRUE
)

# convert DESeq2 results to a regular data frame
# baseMean: mean normalized count across all samples
# log2FoldChange: estimated log2(SAT2 / MOCK)
# pvalue: raw Wald test p value
# padj: BH-adjusted p value across all tested genes
res_sat2_df <- as.data.frame(res_sat2) %>%
  rownames_to_column("gene_name") %>%
  select(gene_name, baseMean, log2FoldChange, pvalue, padj)

# add the per-condition mean normalized counts for easier biological interpretation
fc_table_sat2 <- res_sat2_df %>%
  left_join(group_means, by = "gene_name")

# genes with positive log2FC are upregulated in SAT2 relative to MOCK
fc_sat2_up <- fc_table_sat2 %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange))

# genes with negative log2FC are downregulated in SAT2 relative to MOCK
fc_sat2_down <- fc_table_sat2 %>%
  filter(log2FoldChange < 0) %>%
  arrange(log2FoldChange)

# write all SAT2 comparison outputs
# --- UNCOMMENT THESE LINES TO SAVE CSVs TO WORKING DIRECTORY --- 
# write_csv(fc_table_sat2, "sat2_vs_mock_all_genes.csv")
# write_csv(fc_sat2_up, "sat2_upregulated_genes_with_pvals.csv")
# write_csv(fc_sat2_down, "sat2_downregulated_genes_with_pvals.csv")

# -----------------------------
# M1 vs MOCK
# -----------------------------
# same logic as above, but now comparing M1 to MOCK
# log2FoldChange here is log2(M1 / MOCK)
# pvalue is again from the Wald test on that model coefficient
res_m1 <- results(
  dds,
  contrast             = c("condition", "M1", "MOCK"),
  alpha                = 0.05,
  pAdjustMethod        = "BH",
  independentFiltering = TRUE
)

res_m1_df <- as.data.frame(res_m1) %>%
  rownames_to_column("gene_name") %>%
  select(gene_name, baseMean, log2FoldChange, pvalue, padj)

fc_table_m1 <- res_m1_df %>%
  left_join(group_means, by = "gene_name")

# positive log2FC = up in M1 relative to MOCK
fc_m1_up <- fc_table_m1 %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange))

# negative log2FC = down in M1 relative to MOCK
fc_m1_down <- fc_table_m1 %>%
  filter(log2FoldChange < 0) %>%
  arrange(log2FoldChange)

# write all M1 comparison outputs
# --- UNCOMMENT THESE LINES TO SAVE CSVs TO WORKING DIRECTORY --- 
# write_csv(fc_table_m1, "m1_vs_mock_all_genes.csv")
# write_csv(fc_m1_up, "m1_upregulated_genes_with_pvals.csv")
# write_csv(fc_m1_down, "m1_downregulated_genes_with_pvals.csv")
