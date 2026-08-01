#!/usr/bin/env bash
set -euo pipefail

mkdir -p sequence_processing_output/metaxa2 sequence_processing_output/lotus3_us sequence_processing_output/lotus3_dk

#---(a) Metaxa2 Screening
mapfile -t US_SAMPLES < US_mapped.txt
mapfile -t DK_SAMPLES < DK_mapped.txt

echo "=== Metaxa2 screening: US cohort ==="
: > sequence_processing_output/metaxa2/us_archaea_counts.tsv
for sample_id in "${US_SAMPLES[@]}"; do
    metaxa2 \
        -1 "raw_fastq/us_cohort/${sample_id}_R1.fastq.gz" \
        -2 "raw_fastq/us_cohort/${sample_id}_R2.fastq.gz" \
        -o "sequence_processing_output/metaxa2/${sample_id}" \
        --cpu 4 \
        --reliability 1.0

    n_archaea=$(awk -F'\t' '$2 ~ /^Archaea/' "sequence_processing_output/metaxa2/${sample_id}.taxonomy.txt" | wc -l)
    echo -e "${sample_id}\t${n_archaea}" >> sequence_processing_output/metaxa2/us_archaea_counts.tsv
done

echo "=== Metaxa2 screening: DK cohort ==="
: > sequence_processing_output/metaxa2/dk_archaea_counts.tsv
for sample_id in "${DK_SAMPLES[@]}"; do
    metaxa2 \
        -1 "raw_fastq/dk_cohort/${sample_id}_R1.fastq.gz" \
        -2 "raw_fastq/dk_cohort/${sample_id}_R2.fastq.gz" \
        -o "sequence_processing_output/metaxa2/${sample_id}" \
        --cpu 4 \
        --reliability 1.0

    n_archaea=$(awk -F'\t' '$2 ~ /^Archaea/' "sequence_processing_output/metaxa2/${sample_id}.taxonomy.txt" | wc -l)
    echo -e "${sample_id}\t${n_archaea}" >> sequence_processing_output/metaxa2/dk_archaea_counts.tsv
done

# Combine both cohorts: samples with non-zero archaeal read count -> archaea.txt
awk -F'\t' '$2 > 0 {print $1}' \
    sequence_processing_output/metaxa2/us_archaea_counts.tsv \
    sequence_processing_output/metaxa2/dk_archaea_counts.tsv \
    > archaea.txt

echo "Archaea-positive samples written to archaea.txt ($(wc -l < archaea.txt) samples)"
echo "Expected from paper: 8 DK-positive, 314 US-positive -- verify your counts against this."

#---(b) Lotus3/ DADA2 ASV calling

echo "=== Lotus3 ASV calling: US cohort ==="
cd sequence_processing_output/lotus3_us
lotus3 -i ../../raw_fastq/us_cohort --m lotus_us.txt -p miSeq -amplicon_type SSU \
    -refDB KSGP -taxAligner 2 -o output_lotus3_ASVs \
    clustering dada2 -t 4
cd ../../

echo "=== Lotus3 ASV calling: DK cohort ==="
cd sequence_processing_output/lotus3_dk
lotus3 -i ../../raw_fastq/dk_cohort --m lotus_dk.txt -p miSeq -amplicon_type SSU \
    -refDB KSGP -taxAligner 2 -o output_lotus3_ASVs \
    clustering dada2 -t 4
cd ../../

#---(c) use Lotus3's own family-level table

cp sequence_processing_output/lotus3_us/output_lotus3_ASVs/higherLvl/family.txt \
   us_counts_after_preprocess.tsv

cp sequence_processing_output/lotus3_dk/output_lotus3_ASVs/higherLvl/family.txt \
   de_counts_after_preprocess.tsv

echo "Done. archaea.txt, us_counts_after_preprocess.tsv, and de_counts_after_preprocess.tsv"
echo "Please continue with preprocess.py."
