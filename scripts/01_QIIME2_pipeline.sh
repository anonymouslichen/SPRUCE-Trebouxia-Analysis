#!/bin/bash

################################################################################
# QIIME2 Pipeline for Evernia Photobiont Metabarcoding Analysis
################################################################################
# Project: SPRUCE Trebouxia photobiont community analysis
# Author: Abigail Meyer
# Last edit: 05/12/2022
# 
# Purpose: Process Illumina MiSeq ITS amplicon data to identify and quantify
#          Trebouxia algal photobiont diversity in Evernia mesomorpha lichens
#          from experimental warming treatments
#
# Sequencing: 2x300bp paired-end Illumina MiSeq
# Target: ITS1-ITS2 region using algal-specific primers
# Primers: ITS1T (GGAAGGATCATTGAATCTATCGT) and ITS2T
#
# QIIME2 version: 2021.4
################################################################################

################################################################################
# PART 1: IMPORT RAW SEQUENCING DATA
################################################################################
# Import FASTQ files into QIIME2 format (.qza)
# Using manifest files that specify sample IDs and FASTQ file paths

echo "=== Importing Fall 2020 samples ==="
# Fall 2020 sequencing run - initial collection
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path ~/Evernia_F20/evernia_mappingfile_F20.txt \
  --output-path ~/Evernia_Intermediary_Files/QIIME_files/evernia-f20-single-end-import.qza \
  --input-format SingleEndFastqManifestPhred33V2

echo "=== Importing Spring 2022 samples ==="
# Spring 2022 sequencing run - additional samples
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path ~/Evernia_Folmania_S22/spruce_evernia_mappingfile_S22.txt \
  --output-path ~/Evernia_Intermediary_Files/QIIME_files/evernia-s22-single-end-import.qza \
  --input-format SingleEndFastqManifestPhred33V2

# Move to working directory for subsequent steps
cd ~/Evernia_Intermediary_Files/QIIME_files

################################################################################
# PART 2: REMOVE PRIMER SEQUENCES
################################################################################
# Use cutadapt to remove forward primer (ITS1T) from 5' end of reads
# --p-error-rate 0: No mismatches allowed in primer sequence (stringent)
# This ensures only reads with correct primers are retained

echo "=== Trimming primers from Fall 2020 samples ==="
# R1 primer: GGAAGGATCATTGAATCTATCGT (ITS1T)
qiime cutadapt trim-single \
  --i-demultiplexed-sequences evernia-f20-single-end-import.qza \
  --p-front GGAAGGATCATTGAATCTATCGT \
  --p-error-rate 0 \
  --o-trimmed-sequences evernia-f20-primer-trimmed-single-end.qza \
  --verbose

echo "=== Trimming primers from Spring 2022 samples ==="
qiime cutadapt trim-single \
  --i-demultiplexed-sequences evernia-s22-single-end-import.qza \
  --p-front GGAAGGATCATTGAATCTATCGT \
  --p-error-rate 0 \
  --o-trimmed-sequences evernia-s22-primer-trimmed-single-end.qza \
  --verbose

################################################################################
# PART 3: QUALITY CONTROL VISUALIZATION
################################################################################
# Create quality score visualizations to determine truncation parameters
# View these .qzv files at https://view.qiime2.org/

echo "=== Creating quality visualizations ==="
qiime demux summarize \
  --i-data evernia-f20-primer-trimmed-single-end.qza \
  --o-visualization evernia-f20-primer-trimmed-single-end.qzv

qiime demux summarize \
  --i-data evernia-s22-primer-trimmed-single-end.qza \
  --o-visualization evernia-s22-primer-trimmed-single-end.qzv

################################################################################
# PART 4: DENOISING WITH DADA2
################################################################################
# DADA2 performs:
# 1. Quality filtering
# 2. Denoising (error correction)
# 3. Chimera removal
# 4. Generation of Amplicon Sequence Variants (ASVs)
#
# --p-trunc-len 240: Truncate reads at 240bp (based on quality score inspection)
#                    Chosen where median PHRED score drops below 20

echo "=== Denoising Fall 2020 samples with DADA2 ==="
qiime dada2 denoise-single \
  --i-demultiplexed-seqs evernia-f20-primer-trimmed-single-end.qza \
  --p-trunc-len 240 \
  --o-table evernia_f20_table.qza \
  --o-representative-sequences evernia_f20_rep-seqs.qza \
  --o-denoising-stats evernia_f20_denoising-stats.qza

echo "=== Denoising Spring 2022 samples with DADA2 ==="
qiime dada2 denoise-single \
  --i-demultiplexed-seqs evernia-s22-primer-trimmed-single-end.qza \
  --p-trunc-len 240 \
  --o-table evernia_s22_table.qza \
  --o-representative-sequences evernia_s22_rep-seqs.qza \
  --o-denoising-stats evernia_s22_denoising-stats.qza

################################################################################
# PART 5: MERGE SEQUENCING RUNS
################################################################################
# Combine ASV tables and representative sequences from both sequencing runs
# This creates a single unified dataset for downstream analysis

echo "=== Merging ASV tables from both runs ==="
qiime feature-table merge \
  --i-tables evernia_f20_table.qza \
  --i-tables evernia_s22_table.qza \
  --o-merged-table spruce-evernia-table.qza

echo "=== Merging representative sequences from both runs ==="
qiime feature-table merge-seqs \
  --i-data evernia_f20_rep-seqs.qza \
  --i-data evernia_s22_rep-seqs.qza \
  --o-merged-data spruce-evernia-rep-seqs.qza

################################################################################
# PART 6: VERIFY MERGE SUCCESS
################################################################################
# Create visualizations to confirm successful merging

echo "=== Creating merged data visualizations ==="
qiime feature-table summarize \
  --i-table spruce-evernia-table.qza \
  --o-visualization spruce-evernia-table.qzv

qiime feature-table tabulate-seqs \
  --i-data spruce-evernia-rep-seqs.qza \
  --o-visualization spruce-evernia-rep-seqs.qzv

################################################################################
# PART 7: TAXONOMIC ASSIGNMENT VIA OPEN-REFERENCE CLUSTERING
################################################################################
# Cluster ASVs against a reference database of known Trebouxia sequences
# from Muggia et al. 2020 (comprehensive Trebouxia phylogeny)
#
# Open-reference approach:
# 1. First tries to match to reference sequences at 97% similarity
# 2. Novel sequences that don't match become new reference OTUs
#
# --p-perc-identity 0.97: 97% similarity threshold (standard for species-level)

echo "=== Clustering against Muggia et al. 2020 Trebouxia reference database ==="
qiime vsearch cluster-features-open-reference \
  --i-sequences spruce-evernia-rep-seqs.qza \
  --i-table spruce-evernia-table.qza \
  --i-reference-sequences A_C_I_S_unaligned_ITS.qza \
  --p-perc-identity 0.97 \
  --o-clustered-table trebouxia-table-or-97.qza \
  --o-clustered-sequences trebouxia-rep-seqs-or-97.qza \
  --o-new-reference-sequences trebouxia-new-ref-seqs-or-97.qza

################################################################################
# PART 8: ADD METADATA AND CREATE VISUALIZATIONS
################################################################################
# Incorporate sample metadata (treatment, temperature, etc.) into OTU table

echo "=== Adding metadata to OTU table ==="
qiime metadata tabulate \
  --m-input-file evernia_metadata_S22.txt \
  --m-input-file trebouxia-table-or-97.qza \
  --o-visualization trebouxia-table-or-97-metadata.qzv

echo "=== Creating final OTU table visualizations ==="
qiime feature-table summarize \
  --i-table trebouxia-table-or-97.qza \
  --o-visualization trebouxia-table-or-97.qzv

qiime feature-table tabulate-seqs \
  --i-data trebouxia-rep-seqs-or-97.qza \
  --o-visualization trebouxia-rep-seqs-or-97.qzv

################################################################################
# PART 9: ALPHA RAREFACTION ANALYSIS
################################################################################
# Create rarefaction curves to assess sampling depth adequacy
# --p-max-depth 20000: Maximum rarefaction depth to test

echo "=== Generating alpha rarefaction curves ==="
qiime diversity alpha-rarefaction \
  --i-table trebouxia-table-or-97.qza \
  --p-max-depth 20000 \
  --o-visualization trebouxia-table-or-97-rarefaction.qzv

################################################################################
# PART 10: FILTER LOW-ABUNDANCE OTUs
################################################################################
# Remove rare OTUs that may represent sequencing errors or contaminants
# Threshold: <1047 reads total (0.10% rule from Reitmeier et al. 2021)
#
# This conservative filtering reduces noise while retaining biologically
# relevant diversity

echo "=== Filtering OTUs with <1047 total reads (0.10% threshold) ==="
qiime feature-table filter-features \
  --i-table trebouxia-table-or-97.qza \
  --p-min-frequency 1047 \
  --o-filtered-table trebouxia-table-or-97-low-abundance-remove.qza

################################################################################
# PART 11: FILTER SEQUENCES TO MATCH FILTERED TABLE
################################################################################
# Remove representative sequences for filtered-out OTUs
# Ensures sequence and table files match

echo "=== Removing sequences for filtered OTUs ==="
qiime feature-table filter-seqs \
  --i-data trebouxia-rep-seqs-or-97.qza \
  --i-table trebouxia-table-or-97-low-abundance-remove.qza \
  --o-filtered-data trebouxia-rep-seqs-or-97-low-abundance-remove.qza

################################################################################
# PART 12: FINAL VISUALIZATIONS
################################################################################
# Create visualizations of filtered data

echo "=== Creating filtered data visualizations ==="
qiime feature-table summarize \
  --i-table trebouxia-table-or-97-low-abundance-remove.qza \
  --o-visualization trebouxia-table-or-97-low-abundance-remove.qzv

qiime feature-table tabulate-seqs \
  --i-data trebouxia-rep-seqs-or-97-low-abundance-remove.qza \
  --o-visualization trebouxia-rep-seqs-or-97-low-abundance-remove.qzv

################################################################################
# PART 13: EXPORT FOR DOWNSTREAM ANALYSIS IN R
################################################################################
# Export OTU table from QIIME2 format (.qza) to standard formats
# BIOM format is then converted to tab-separated text for R

echo "=== Exporting OTU table ==="
qiime tools export \
  --input-path trebouxia-table-or-97-low-abundance-remove.qza \
  --output-path ~/Evernia_Intermediary_Files/QIIME_files

# Convert BIOM format to tab-separated values (TSV)
echo "=== Converting BIOM to TSV format ==="
biom convert \
  -i feature-table.biom \
  -o spruce-evernia-feature-table.tsv \
  --to-tsv

echo "=== QIIME2 pipeline complete! ==="
echo "Output files ready for R analysis"

################################################################################
# KEY OUTPUT FILES:
# - trebouxia-table-or-97-low-abundance-remove.qza: Final OTU table
# - trebouxia-rep-seqs-or-97-low-abundance-remove.qza: Representative sequences
# - spruce-evernia-feature-table.tsv: OTU table in text format for R
# - Multiple .qzv visualization files (view at https://view.qiime2.org/)
################################################################################
