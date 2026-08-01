# permanova_maaslin2.R
# -----------------------------------------------------------------------------
# Implements Methods Sections 3.4 (batch-effect characterization), 4.1 (primary
# GLM), 4.2 (PERMANOVA-based beta-diversity), and 4.3 (MaAsLin2 mixed-effects
# model). Reads the outputs of preprocess.py (permanova/R_meta.tsv,
# permanova/concated_clr.tsv) and writes archaea_sig.tsv / preterm_sig.tsv in
# the format already consumed by analysis.py.
#
# PROVENANCE NOTE: this script reconstructs the described methodology from the
# manuscript text. It has not been run against the actual cohort data, so it
# is not verified to reproduce the exact reported statistics (beta = 0.457,
# R^2 ~= 0.01364, pseudo-F = 31.089, etc.). Run it against your real data and
# confirm those numbers match before treating this as the authoritative
# analysis code.
#
# Dependencies: vegan, Maaslin2 (Bioconductor), tidyverse (readr/dplyr for
# convenience; not strictly required)
# -----------------------------------------------------------------------------

library(vegan)
library(Maaslin2)
library(dplyr)
library(readr)

dir.create("permanova", showWarnings = FALSE)

# ----------------------------- 1) Load inputs --------------------------------
meta <- read_tsv("permanova/R_meta.tsv") %>% rename(sample_id = 1)
clr  <- read_tsv("permanova/concated_clr.tsv") %>% rename(sample_id = 1)

stopifnot(all(meta$sample_id == clr$sample_id))

taxa_cols <- setdiff(colnames(clr), "sample_id")
clr_mat <- as.matrix(clr[, taxa_cols])
rownames(clr_mat) <- clr$sample_id

meta$Archaea <- as.integer(meta$Archaea)
meta$Preterm <- as.integer(meta$Preterm)
meta$cohort  <- factor(meta$cohort, levels = c("US", "DK"))
meta$PatID   <- factor(meta$PatID)

# ----------------------------- 2) Primary GLM (Section 4.1) ------------------
# Delivery mode ~ archaeal presence, adjusting for cohort.
glm_fit <- glm(Preterm ~ Archaea + cohort, data = meta, family = binomial())
glm_summary <- summary(glm_fit)
cat("=== Primary GLM: Preterm ~ Archaea + cohort ===\n")
print(glm_summary$coefficients)
# Expect (per manuscript): beta (Archaea) = 0.457, SE = 0.130, p < 0.001
write.csv(glm_summary$coefficients, "primary_glm_results.csv")

# ----------------------------- 3) Batch-effect PERMANOVA (Section 3.4) -------
dist_mat <- vegdist(clr_mat, method = "euclidean")  # Euclidean on CLR == Aitchison

cat("=== Batch-effect PERMANOVA (cohort) ===\n")
adonis_cohort <- adonis2(dist_mat ~ cohort, data = meta, permutations = 999)
print(adonis_cohort)
# Expect (per manuscript): R^2 ~= 0.013, pseudo-F = 28.472, p = 0.001

permdisp_cohort <- betadisper(dist_mat, meta$cohort)
print(permutest(permdisp_cohort, permutations = 999))

# ----------------------------- 4) Archaea / birth-outcome PERMANOVA (4.2) ----
# Permutations restricted within subjects (strata = PatID) to respect the
# longitudinal/repeated-measures structure.
cat("=== PERMANOVA: Archaea presence ===\n")
adonis_archaea <- adonis2(
  dist_mat ~ Archaea, data = meta, permutations = 999, strata = meta$PatID
)
print(adonis_archaea)
# Expect: R^2 ~= 0.01364, pseudo-F = 31.089, p = 0.001

cat("=== PERMANOVA: Birth outcome ===\n")
adonis_preterm <- adonis2(
  dist_mat ~ Preterm, data = meta, permutations = 999, strata = meta$PatID
)
print(adonis_preterm)
# Expect: R^2 ~= 0.0104, pseudo-F = 23.616, p = 0.001

# ----------------------------- 5) PERMDISP (Fig 3) ---------------------------
disp_archaea <- betadisper(dist_mat, factor(meta$Archaea, labels = c("Absent", "Present")))
disp_preterm <- betadisper(dist_mat, factor(meta$Preterm, labels = c("Term", "Preterm")))

dist_to_centroid_df <- data.frame(
  sample_id = names(disp_archaea$distances),
  archaea_group = disp_archaea$group,
  preterm_group = disp_preterm$group,
  distance_to_centroid = disp_archaea$distances
)
write_tsv(dist_to_centroid_df, "permanova/distance_to_centroid.tsv")

wilcox_archaea <- wilcox.test(distances ~ group, data = data.frame(
  distances = disp_archaea$distances, group = disp_archaea$group
))
wilcox_preterm <- wilcox.test(distances ~ group, data = data.frame(
  distances = disp_preterm$distances, group = disp_preterm$group
))
cat("PERMDISP Wilcoxon (archaea):\n"); print(wilcox_archaea)
cat("PERMDISP Wilcoxon (preterm):\n"); print(wilcox_preterm)

# ----------------------------- 6) MaAsLin2 mixed-effects model (4.3) --------
# Eq. 7-9: y_itj = b0 + b1*Preterm + b2*cohort + b3*Archaea + b_i (random) + e
maaslin_input_data <- as.data.frame(clr_mat)
maaslin_input_meta <- as.data.frame(meta)
rownames(maaslin_input_meta) <- meta$sample_id

fit <- Maaslin2(
  input_data      = maaslin_input_data,
  input_metadata  = maaslin_input_meta,
  output          = "maaslin2_output",
  fixed_effects   = c("Preterm", "cohort", "Archaea"),
  random_effects  = c("PatID"),
  normalization   = "NONE",   # data is already CLR-transformed
  transform       = "NONE",
  standardize     = FALSE,
  min_prevalence  = 0,
  correction      = "BH"
)

all_results <- fit$results

# ----------------------------- 7) Format sig tables to match existing files -
# Reproduces the exact column layout already used by archaea_sig.tsv /
# preterm_sig.tsv: feature, metadata, direction, N, N_not_zero, pval, qval.
# "direction" is derived from the sign of the coefficient; verify the
# reference level for each factor (Archaea: 0 = absent is reference; Preterm:
# 0 = term is reference) matches how your metadata columns are actually coded
# before trusting the direction labels below.
format_sig_table <- function(results, metadata_name, positive_label, negative_label, fdr_threshold = 0.05) {
  results %>%
    filter(metadata == metadata_name, qval < fdr_threshold) %>%
    mutate(direction = ifelse(coef > 0, positive_label, negative_label)) %>%
    select(feature, metadata, direction, N, N.not.0, pval, qval) %>%
    rename(N_not_zero = N.not.0) %>%
    arrange(qval)
}

archaea_sig <- format_sig_table(all_results, "Archaea", "present", "absent")
preterm_sig <- format_sig_table(all_results, "Preterm", "preterm", "term")

write_csv(archaea_sig, "archaea_sig.tsv")
write_csv(preterm_sig, "preterm_sig.tsv")

cat(sprintf(
  "Wrote archaea_sig.tsv (%d families, FDR<0.05) and preterm_sig.tsv (%d families, FDR<0.05)\n",
  nrow(archaea_sig), nrow(preterm_sig)
))
# Expect (per manuscript): 42 families for Archaea, 13 families for Preterm
