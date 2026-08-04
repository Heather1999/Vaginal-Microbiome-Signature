# Vaginal-Microbiome-Signature

Analysis code accompanying "Possible association of Vaginal archaeal presence
with preterm birth and a distinct anaerobic microbiome configuration."
Integrates a US discovery cohort (Callahan et al. 2017, n=2179) and a Danish
validation cohort (Borum et al. 2025, n=71) to test whether vaginal archaeal
presence is associated with preterm birth.

## Pipeline order

```
0. (not scripted)                 Cohort clinical/participant metadata --
                                   us_meta_after_preprocess.tsv and 
                                   de_meta_after_preprocess.tsv -- 
                                   sourced from the original Callahan et 
                                   al. 2017 public dataset and the Borum 
                                   et al. 2025 trial data respectively 
                                   and processed from scripts that are 
                                   not available here.

1. sequence_processing.sh         Metaxa2 archaeal screening (Sec. 2) +
                                   Lotus3/DADA2 ASV calling (Sec. 3.1).
                                   Reads sample IDs from US_mapped.txt /
                                   DK_mapped.txt. Uses Lotus3's own
                                   higherLvl/family.txt output directly --
                                   no separate aggregation script needed.
                                   outputs: archaea.txt,
                                            us_counts_after_preprocess.tsv,
                                            de_counts_after_preprocess.tsv

2. preprocess.py                  CLR transformation + harmonized US/DK
                                   metadata (Sec. 3.2-3.3)
                                   outputs: permanova/R_meta.tsv,
                                            permanova/concated_clr.tsv,
                                            permanova/dist.tsv,
                                            permanova/pca_us_de.pdf (Fig 1)

3. permanova_maaslin2.R           Batch-effect + archaea/birth-outcome
                                   PERMANOVA and PERMDISP (Sec. 3.4, 4.2,
                                   Fig 3), primary GLM (Sec. 4.1), MaAsLin2
                                   mixed-effects model (Sec. 4.3)
                                   outputs: archaea_sig.tsv, preterm_sig.tsv,
                                            primary_glm_results.csv

4. prediction_archaea.py          Elastic net logistic regression, trained
   prediction_birth.py            on US / externally validated on DK
                                   (Sec. 5, Fig 4)
                                   outputs: Strict_US_trained_Features_Archaea.tsv,
                                            Strict_US_trained_Features.tsv,
                                            archaea_model_report.pdf,
                                            birth_model_report.pdf

5. analysis.py                    Overlap between MaAsLin2- and elastic-net-
                                   selected taxa (Fig 5)
```

Section 6 of the manuscript (shotgun metagenomic confirmation on the 8
Danish archaea-positive samples: NGLess preprocessing, Kraken2 classification,
targeted Methanobrevibacter genome mapping against GTDB R226) is not yet
represented in this repository.

## Dependencies

Python: see `requirements.txt` (`pip install -r requirements.txt`).

R: see `install_r_dependencies.R` (`vegan`, `Maaslin2` via Bioconductor,
`dplyr`, `readr`).

External tools (not Python/R packages, must be installed separately):
Metaxa2 2.2.3, Lotus3 3.2.2.

## Provenance / verification status

`sequence_processing.sh`, `permanova_maaslin2.R`, and the bug fixes in
`preprocess.py` were reconstructed from the manuscript's Methods section to
fill gaps in this repository. They have not been run against the original
raw data, so they are not yet verified to reproduce the paper's exact
reported numbers. Before treating them as final:

- Confirm `archaea.txt` yields 314 US-positive and 8 DK-positive samples
  (Section 2.1).
- Confirm the batch-effect PERMANOVA gives R^2 ~= 0.013, pseudo-F = 28.472,
  p = 0.001 (Section 3.4).
- Confirm the archaea/birth-outcome PERMANOVA give R^2 ~= 0.01364 /
  0.0104 and pseudo-F = 31.089 / 23.616 respectively (Section 4.2).
- Confirm the primary GLM gives beta = 0.457, SE = 0.130, p < 0.001
  (Section 4.1).
- Confirm MaAsLin2 yields 42 (Archaea) and 13 (Preterm) FDR<0.05 families
  (Section 4.3) matching the existing `archaea_sig.tsv` / `preterm_sig.tsv`.

If any of these don't match, the most likely culprits are the factor
reference-level assumptions in `permanova_maaslin2.R`'s `direction` column,
or the taxonomy rank naming conventions in Lotus3's `higherLvl/family.txt`
output for your installed KSGP database version.
