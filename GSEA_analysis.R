library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(DT)

setwd("C:/Users/lisam/Documents/RNA_Seq_2026/gene_expression")

# -------------------------------------------------------
# 1. BUILD RANKED GENE LIST from DESeq2 SAT2 result
# -------------------------------------------------------
# fc_table_sat2 already exists from  DGE script
# rank by sign(log2FC) * -log10(pvalue) to capture
# direction and significance
# NOTE: replace any pvalue = 0 with smallest nonzero double
# to avoid Inf in rank metric (common with DESeq2)

ranked_df <- fc_table_sat2 %>%
  filter(!is.na(pvalue), !is.na(log2FoldChange)) %>%
  mutate(
    pvalue_safe  = ifelse(pvalue == 0, .Machine$double.xmin, pvalue),
    rank_metric  = sign(log2FoldChange) * -log10(pvalue_safe)
  )

mydata.gsea <- ranked_df$rank_metric
names(mydata.gsea) <- ranked_df$gene_name
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

# look at top and bottom of ranked list
head(mydata.gsea, 10)
tail(mydata.gsea, 10)

# -------------------------------------------------------
# 2. BUILD GENE SETS
# -------------------------------------------------------

# --- Hallmark: EMT, IFN-gamma, IFN-alpha, DNA repair, heat shock, etc. ---
hs_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

# --- C2 curated: PRC2 targets, interferon, TGFb, DNA damage ---
hs_c2 <- msigdbr(species = "Homo sapiens", category = "C2") %>%
  dplyr::select(gs_name, gene_symbol)

# --- Custom HSAT2 gene set from literature ---
# Ponomartsev et al. 2023 (IJMS): EMT related genes
hsat2_custom <- tibble(
  gs_name = c(
    rep("HSAT2_EMT_DIRECT_UP", 7),
    rep("HSAT2_EMT_DIRECT_DOWN", 1)
  ),
  gene_symbol = c(
    "SNAI1", "SNAI2", "ZEB1", "VIM", "ACTA2", "COL1A1", "COL11A1",  # upregulated by HSAT2
    "CDH1"                                                          # downregulated by HSAT2
  )
)

# combine all sets into a tibble
all_sets <- bind_rows(hs_hallmark, hs_c2, hsat2_custom)

# -------------------------------------------------------
# 3. RUN GSEA — SAT2 vs MOCK (clusterProfiler)
# -------------------------------------------------------
# pvalueCutoff = 1 keeps all results so we can filter ourselves
# NES > 0 = pathway enriched when HSAT2 is knocked down
# NES < 0 = pathway depleted when HSAT2 is knocked down

set.seed(123)
myGSEA.res <- GSEA(
  mydata.gsea,
  TERM2GENE    = all_sets,
  pvalueCutoff = 1,
  verbose      = FALSE
)

myGSEA.df <- as_tibble(myGSEA.res@result) %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "suppressed_by_HSAT2",  # upregulated after KD, so HSAT2 suppresses it
    NES < 0 ~ "sustained_by_HSAT2"    # downregulated after KD, so HSAT2 sustains it
  ))

# -------------------------------------------------------
# 4. VIEW RESULTS — interactive searchable table
# -------------------------------------------------------
datatable(myGSEA.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption    = 'GSEA results: SAT2 gapmer vs MOCK',
          options    = list(keys = TRUE, searchHighlight = TRUE,
                            pageLength = 10,
                            lengthMenu = c("10", "25", "50"))) %>%
  formatRound(columns = c(2:10), digits = 2)

# -------------------------------------------------------
# 5. LOOK AT HSAT2-RELEVANT PATHWAYS
# -------------------------------------------------------
myGSEA.df %>%
  filter(str_detect(ID,
                    "EMT|INTERFERON|PRC2|SUZ12|HSAT2|TGF|HSP|DNA_REPAIR|HOMOLOGOUS")) %>%
  dplyr::select(ID, NES, p.adjust, setSize, phenotype) %>%
  arrange(p.adjust) %>%
  print(n = 30)

# -------------------------------------------------------
# 6. CHECK POSITIVE NES — what does HSAT2 normally suppress?
# -------------------------------------------------------
myGSEA.df %>%
  filter(NES > 0, p.adjust < 0.05) %>%
  arrange(p.adjust) %>%
  dplyr::select(ID, NES, p.adjust, setSize) %>%
  print(n = 30)

# -------------------------------------------------------
# 7. RUN GSEA — M1 vs MOCK (positive control)
# -------------------------------------------------------
# Any pathway significant in BOTH SAT2 and M1 is likely not HSAT2-specific 

ranked_m1 <- fc_table_m1 %>%
  filter(!is.na(pvalue), !is.na(log2FoldChange)) %>%
  mutate(
    pvalue_safe = ifelse(pvalue == 0, .Machine$double.xmin, pvalue),
    rank_metric = sign(log2FoldChange) * -log10(pvalue_safe)
  )

m1.gsea <- setNames(ranked_m1$rank_metric, ranked_m1$gene_name)
m1.gsea <- sort(m1.gsea, decreasing = TRUE)

set.seed(123)
m1GSEA.res <- GSEA(
  m1.gsea,
  TERM2GENE    = all_sets,
  pvalueCutoff = 1,
  verbose      = FALSE
)

m1GSEA.df <- as_tibble(m1GSEA.res@result) %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "suppressed_by_M1",
    NES < 0 ~ "sustained_by_M1"
  ))

# -------------------------------------------------------
# 8. COMPARE SAT2 vs M1 — identify HSAT2-specific hits
# -------------------------------------------------------
sat2_sig_ids <- myGSEA.df %>% filter(p.adjust < 0.05) %>% pull(ID)
m1_sig_ids   <- m1GSEA.df  %>% filter(p.adjust < 0.05) %>% pull(ID)

# SAT2-specific pathways (not significant in M1 control) — highest confidence
hsat2_specific <- setdiff(sat2_sig_ids, m1_sig_ids)
cat("Pathways specific to HSAT2 knockdown (not in M1 control):\n")
print(hsat2_specific)

# Full side-by-side NES comparison for all shared significant pathways
sat2_nes <- myGSEA.df %>% filter(p.adjust < 0.05) %>% dplyr::select(ID, NES, p.adjust)
m1_nes   <- m1GSEA.df  %>% filter(p.adjust < 0.05) %>% dplyr::select(ID, NES, p.adjust)

comparison <- left_join(sat2_nes, m1_nes, by = "ID", suffix = c("_SAT2", "_M1"))
print(comparison)

# pathways where SAT2 and M1 go in OPPOSITE directions are especially interesting
comparison %>%
  filter(!is.na(NES_M1)) %>%
  filter(sign(NES_SAT2) != sign(NES_M1)) %>%
  arrange(p.adjust_SAT2)

myGSEA.df %>%
  filter(p.adjust < 0.05) %>%
  # write_csv("gsea_sat2_vs_mock_significant.csv")

m1GSEA.df %>%
  filter(p.adjust < 0.05) %>%
  # write_csv("gsea_m1_vs_mock_significant.csv")

# positive NES for SAT2 — pathways normally SUPPRESSED by HSAT2
sat2_positive <- myGSEA.df %>%
  filter(NES > 0, p.adjust < 0.05) %>%
  arrange(desc(NES))

# write_csv(sat2_positive, "gsea_sat2_positive_nes_significant.csv")

# positive NES for M1
m1_positive <- m1GSEA.df %>%
  filter(NES > 0, p.adjust < 0.05) %>%
  arrange(desc(NES))

# write_csv(m1_positive, "gsea_m1_positive_nes_significant.csv")