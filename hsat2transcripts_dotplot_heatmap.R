library(tidyverse)
library(tximport)
library(edgeR)
library(cowplot)
library(matrixStats)
library(pheatmap)

# work in gene_expression folder
setwd("C:/Users/lisam/Documents/RNA_Seq_2026/gene_expression")

# where the kallisto output folders live
kallisto_dir <- "C:/Users/lisam/Documents/RNA_Seq_2026/kallisto_outputs_t2t_48hr_gapmer"

# read study design and group columns by treatment
targets <- read_tsv("studydesign.txt", col_types = cols()) %>%
  mutate(
    Treatment = factor(toupper(Treatment), levels = c("MOCK", "M1", "SAT2")),
    Replicates = factor(as.character(Replicates), levels = c("1", "2", "3", "4"))
  ) %>%
  arrange(Treatment, Replicates)

sample_names <- targets$Name

# point to abundance files
h5_paths  <- file.path(kallisto_dir, sample_names, "abundance.h5")
tsv_paths <- file.path(kallisto_dir, sample_names, "abundance.tsv")
files <- ifelse(file.exists(h5_paths), h5_paths, tsv_paths)
names(files) <- sample_names

# check files exist
missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  stop("Missing abundance files:\n", paste(missing_files, collapse = "\n"))
}

# transcript-level import
txi <- tximport(
  files,
  type = "kallisto",
  txOut = TRUE,
  countsFromAbundance = "lengthScaledTPM",
  ignoreTxVersion = TRUE
)

# exact target HSAT2 transcripts
hsat2_targets <- c(
  "hub_3267197_GCA_009914755.4_hub_3267197_censat_hsat2_17_2(B)",
  "hub_3267197_GCA_009914755.4_hub_3267197_censat_hsat2_16_15(B)",
  "hub_3267197_GCA_009914755.4_hub_3267197_censat_hsat2_10_3(A1,A2)",
  "hub_3267197_GCA_009914755.4_hub_3267197_censat_hsat2_9_1(B)",
  "hub_3267197_GCA_009914755.4_hub_3267197_censat_hsat2_7_6(B)",
  "hub_3267197_GCA_009914755.4_hub_3267197_censat_hsat2_7_12(B)",
  "hub_3267197_GCA_009914755.4_hub_3267197_censat_hsat2_7_13(B)",
  "hub_3267197_GCA_009914755.4_hub_3267197_censat_hsat2_1_3(A1,B)",
  "hub_3267197_GCA_009914755.4_hub_3267197_censat_hsat2_1_4(B,A1,A2)",
  "hub_3267197_GCA_009914755.4_hub_3267197_censat_hsat2_1_6(A2,A1,B)"
)

# -----------------------------
# UNFILTERED matrix for 10-target export
# -----------------------------
dge_all <- DGEList(counts = txi$counts)
dge_all <- calcNormFactors(dge_all, method = "TMM")
cpm.log.all <- cpm(dge_all, log = TRUE)

desired_cols <- targets$Name

full_tbl_all <- as_tibble(cpm.log.all, rownames = "geneID") %>%
  select(geneID, all_of(desired_cols))

hsat2_target_tbl <- full_tbl_all %>%
  filter(geneID %in% hsat2_targets) %>%
  mutate(geneID = factor(geneID, levels = hsat2_targets)) %>%
  arrange(geneID) %>%
  mutate(geneID = as.character(geneID)) %>%
  mutate(
    short_name = str_extract(geneID, "hsat2.*$"),
    MOCK_avg = rowMeans(select(., starts_with("MOCK")), na.rm = TRUE),
    M1_avg   = rowMeans(select(., starts_with("M1")), na.rm = TRUE),
    SAT2_avg = rowMeans(select(., starts_with("SAT2")), na.rm = TRUE),
    M1_vs_MOCK   = M1_avg - MOCK_avg,
    SAT2_vs_MOCK = SAT2_avg - MOCK_avg
  )

write_csv(hsat2_target_tbl, "hsat2_target_transcript_norm_log2cpm_ALL10.csv")

# -----------------------------
# FILTERED matrix for all HSAT2 export
# -----------------------------
keep <- rowSums(cpm(dge_all, log = FALSE) > 1) >= 3
dge_f <- dge_all[keep, , keep.lib.sizes = FALSE]
dge_f <- calcNormFactors(dge_f, method = "TMM")
cpm.log.fn <- cpm(dge_f, log = TRUE)

full_tbl_filtered <- as_tibble(cpm.log.fn, rownames = "geneID") %>%
  select(geneID, all_of(desired_cols))

hsat2_tbl <- full_tbl_filtered %>%
  filter(str_detect(geneID, "hsat2")) %>%
  mutate(short_name = str_extract(geneID, "hsat2.*$"))

write_csv(hsat2_tbl, "hsat2_transcript_filtered_norm_log2cpm.csv")

# -----------------------------
# DOT PLOT: change relative to mock
# -----------------------------
plot_df <- hsat2_target_tbl %>%
  select(short_name, M1_vs_MOCK, SAT2_vs_MOCK) %>%
  pivot_longer(
    cols = c(M1_vs_MOCK, SAT2_vs_MOCK),
    names_to = "Comparison",
    values_to = "log2CPM_diff"
  ) %>%
  mutate(
    Comparison = recode(
      Comparison,
      M1_vs_MOCK = "M1 vs Mock",
      SAT2_vs_MOCK = "SAT2 vs Mock"
    ),
    short_name = factor(short_name, levels = rev(hsat2_target_tbl$short_name))
  )

p_dot <- ggplot(plot_df, aes(x = log2CPM_diff, y = short_name, color = Comparison)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  labs(
    title = "HSAT2 transcript changes relative to mock",
    x = "Difference in mean log2CPM relative to mock",
    y = "HSAT2 transcript"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("hsat2_dotplot_vs_mock.png", p_dot, width = 8, height = 5, dpi = 300)
ggsave("hsat2_dotplot_vs_mock.pdf", p_dot, width = 8, height = 5)

# -----------------------------
# HEATMAP CODE
# -----------------------------
heatmap_mat <- hsat2_target_tbl %>%
  select(short_name, all_of(desired_cols)) %>%
  column_to_rownames("short_name") %>%
  as.matrix()

annotation_col <- targets %>%
  select(Name, Treatment) %>%
  as.data.frame()
rownames(annotation_col) <- annotation_col$Name
annotation_col$Name <- NULL

pheatmap(
  heatmap_mat,
  scale = "none",
  annotation_col = annotation_col,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  fontsize_row = 10,
  main = "HSAT2 target transcript expression across samples",
  filename = "hsat2_heatmap_log2cpm_none_scaling.jpg",
  width = 10,
  height = 7
)

pheatmap(
  heatmap_mat,
  scale = "row",
  annotation_col = annotation_col,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  fontsize_row = 10,
  main = "HSAT2 target transcript expression across samples",
  filename = "hsat2_heatmap_log2cpm_row_scaling.jpg",
  width = 10,
  height = 7
)

# -----------------------------
# DOTPLOT CODE
# -----------------------------

plot_df <- hsat2_target_tbl %>%
  select(short_name, SAT2_vs_MOCK) %>%
  mutate(
    short_name = factor(short_name, levels = rev(hsat2_target_tbl$short_name))
  )

p_dot <- ggplot(plot_df, aes(x = SAT2_vs_MOCK, y = short_name)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3) +
  labs(
    title = "HSAT2 transcript changes after HSAT2 GapmeR treatment",
    x = "Mean log2CPM difference (SAT2 - Mock)",
    y = "HSAT2 transcript"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(hjust = 0.5)
  )

p_dot

ggsave("hsat2_dotplot_sat2_vs_mock.jpg", p_dot, width = 7, height = 5, dpi = 300)

